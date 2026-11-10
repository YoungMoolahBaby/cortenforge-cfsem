//! `steer`: the swimmer `shape` chose, as a rigid body in CortenForge's sim-core: gravity,
//! buoyancy, the magnet's torque and pull, and contact with the jar, with the fluid's forces
//! from our Stokes model through `xfrc_applied` (sim-core's own fluid model can't swim a helix,
//! `check`).
//!
//! The fluid's force is -R [v; w] (R the 6 x 6 resistance matrix about the centre of mass, body
//! frame). The swimmer settles to a new speed in under a millisecond, so a plain explicit force
//! would need steps of microseconds. Instead each step applies the drag of the velocity the step
//! will end with, ignoring contact (linearly implicit): v' = (M + dt R)^-1 (M v + dt F), with F
//! the other known forces (gravity, buoyancy, the magnet). Free of contact this is exact for
//! sim-core's update (v' = v + dt M^-1 (F - R v')), stable at any step and, at steps much longer
//! than the settling time, the overdamped (Stokes) motion. In contact it lags: the drag assumed
//! the swimmer moves freely, so the contact force is too small and friction too weak, unless the
//! step is short (`steer wall`). Section 1 checks the coupling against the model.
//!
//! Section 2: swimming level (the nose droops and swings) and how fast it can pitch and turn, at
//! 2-5 mT. Section 3: a printed spiral staircase in the jar (a post, two fins with windows),
//! stood up off the floor with a still field, then steered by pure pursuit from a camera, its
//! aim corrected for sinking; the field from the coils. `steer course` runs section 3 alone;
//! `steer wall` and `steer contact` probe sim-core's contacts (capsule-box and some cylinder
//! contacts are broken in 0.9, so the swimmer's contact shape is spheres).

use std::f64::consts::PI;
use std::time::Instant;

use cortenforge::sim::core::{Data, Model};
use cortenforge::sim::mjcf::load_model;
use nalgebra::{Matrix3, Matrix6, Quaternion, Rotation3, UnitQuaternion, Vector6};
use rayon::prelude::*;

use crate::check::{GLYCEROL_MU, GLYCEROL_RHO};
use crate::coils::{self, LIQ_H, LIQ_R};
use crate::helix::V3;
use crate::rft::{swim, swim_locked};
use crate::shape;
use crate::swimmer::Swimmer;

const G: f64 = 9.81;
/// Run the field at this fraction of step-out (as in `coils`).
const KAPPA: f64 = 0.75;

fn skew(d: &V3) -> Matrix3<f64> {
    Matrix3::new(0.0, -d.z, d.y, d.z, 0.0, -d.x, -d.y, d.x, 0.0)
}

/// The swimmer as a rigid body. Body frame: origin at the centre of mass, z along the axis
/// (nose up), x along the magnet's moment.
pub struct Body {
    pub mass: f64,
    pub volume: f64,
    /// The centre of buoyancy and the magnet's centre, from the centre of mass.
    pub cob: V3,
    pub mag: V3,
    /// About the centre of mass, body frame.
    pub inertia: Matrix3<f64>,
    pub r: Matrix6<f64>,
    /// (A m^2)
    pub m: f64,
    /// Contact shapes: capsules (from, to, radius), from the centre of mass.
    pub capsules: Vec<(V3, V3, f64)>,
    /// The middle of the body along the axis, from the centre of mass (the point steered).
    pub mid: f64,
    /// The contact shape: spheres (centre from the centre of mass, radius).
    pub spheres: Vec<(V3, f64)>,
    /// The model's predictions, axis held: m/s per Hz, step-out Hz per tesla, sinking upright
    /// and level (m/s); axis free: m/s per Hz, step-out Hz per tesla.
    pub v_hz: f64,
    pub fso_t: f64,
    pub sink_up: f64,
    pub sink_level: f64,
    pub v_hz_free: f64,
    pub fso_t_free: f64,
    pub model_points: usize,
    pub model_secs: f64,
    /// The sampled mass, against the formula's (a check on the inertia's samples).
    pub mass_sampled: f64,
}

impl Body {
    pub fn new(sw: Swimmer) -> Body {
        let ms = sw.mass();
        let com = ms.com;
        let model = sw.model(GLYCEROL_MU);
        // Move the reference point from the swimmer's centre to the centre of mass:
        // [U_c; W] = P [U_p; W], P = [[I, [d]x], [0, I]], d = p - c; R_p = P^T R_c P.
        let d = com - sw.centre();
        let mut p = Matrix6::identity();
        p.fixed_view_mut::<3, 3>(0, 3).copy_from(&skew(&d));
        let r = p.transpose() * model.r * p;
        let held = swim_locked(&model.r);
        let free = swim(&model.r);
        let m = sw.moment();
        let w = (ms.mass - GLYCEROL_RHO * ms.volume) * G;
        let a_inv = r
            .fixed_view::<3, 3>(0, 0)
            .into_owned()
            .try_inverse()
            .expect("drag invertible");
        let (inertia, mass_sampled) = inertia(&sw, &com);
        let caps = capsules(&sw, &com);
        let lo = caps.iter().map(|c| c.0.z.min(c.1.z) - c.2).fold(f64::INFINITY, f64::min);
        let hi = caps.iter().map(|c| c.0.z.max(c.1.z) + c.2).fold(f64::NEG_INFINITY, f64::max);
        let sph = spheres(&caps);
        Body {
            mass: ms.mass,
            volume: ms.volume,
            cob: ms.cob - com,
            mag: sw.magnet_centre() - com,
            inertia,
            r,
            m,
            capsules: caps,
            mid: 0.5 * (lo + hi),
            spheres: sph,
            v_hz: held.u.z * 2.0 * PI,
            fso_t: m / (2.0 * PI * held.torque),
            sink_up: -(a_inv * V3::new(0.0, 0.0, -w)).z,
            sink_level: -(a_inv * V3::new(-w, 0.0, 0.0)).x,
            v_hz_free: free.u.z * 2.0 * PI,
            fso_t_free: m / (2.0 * PI * free.torque),
            model_points: model.points,
            model_secs: model.seconds,
            mass_sampled,
        }
    }

    /// Net weight in glycerol (N).
    pub fn weight(&self) -> f64 {
        (self.mass - GLYCEROL_RHO * self.volume) * G
    }

    /// Tail and nose along the axis, from the centre of mass.
    pub fn ends(&self) -> (f64, f64) {
        let c = &self.capsules;
        (
            c.iter().map(|c| c.0.z.min(c.1.z) - c.2).fold(f64::INFINITY, f64::min),
            c.iter().map(|c| c.0.z.max(c.1.z) + c.2).fold(f64::NEG_INFINITY, f64::max),
        )
    }
}

/// The inertia about the centre of mass from samples: the wire as a line of mass, the head on a
/// fine grid (resin, the magnet denser). Inertia hardly matters here (the swimmer settles in
/// under a millisecond); mass and centre of mass come from `Swimmer::mass`. Returns the inertia
/// and the sampled mass.
fn inertia(sw: &Swimmer, com: &V3) -> (Matrix3<f64>, f64) {
    let h = &sw.helix;
    let mut i = Matrix3::zeros();
    let mut mass = 0.0;
    let mut add = |p: V3, dm: f64| {
        let q = p - com;
        i += (Matrix3::identity() * q.dot(&q) - q * q.transpose()) * dm;
        mass += dm;
    };
    let n = 4000;
    let ds = h.length() / n as f64;
    let lin = sw.rho_body * PI * h.wire * h.wire;
    for k in 0..n {
        add(h.point(h.t_max() * (k as f64 + 0.5) / n as f64), lin * ds);
    }
    let (r, l) = (sw.head_radius, sw.head_len);
    let zc = sw.head_z0() + r;
    let (z0, z1) = (sw.head_z0(), sw.nose_z());
    let step = r / 40.0;
    let mc = sw.magnet_centre();
    let (nx, nz) = ((2.0 * r / step).ceil() as usize, ((z1 - z0) / step).ceil() as usize);
    for a in 0..nx {
        for b in 0..nx {
            for c in 0..nz {
                let p = V3::new(
                    -r + (a as f64 + 0.5) * step,
                    -r + (b as f64 + 0.5) * step,
                    z0 + (c as f64 + 0.5) * step,
                );
                let rho2 = p.x * p.x + p.y * p.y;
                let inside = if p.z < zc {
                    rho2 + (p.z - zc).powi(2) < r * r
                } else if p.z <= zc + l {
                    rho2 < r * r
                } else {
                    !sw.flat_nose && rho2 + (p.z - zc - l).powi(2) < r * r
                };
                if !inside || p.z > z1 {
                    continue;
                }
                let in_mag = rho2 < sw.mag_radius.powi(2) && (p.z - mc.z).abs() < 0.5 * sw.mag_len;
                let rho = if in_mag { sw.rho_magnet } else { sw.rho_body };
                add(p, rho * step.powi(3));
            }
        }
    }
    (i, mass)
}

/// The body's shape for clearances: 16 capsules a turn along the wire, the strut (two), and the
/// head as one capsule reaching the flat nose (its rim rounded, 0.4 mm inside the real corner).
fn capsules(sw: &Swimmer, com: &V3) -> Vec<(V3, V3, f64)> {
    let h = &sw.helix;
    let mut out = Vec::new();
    let n = (16.0 * h.turns) as usize;
    for k in 0..n {
        let a = h.point(h.t_max() * k as f64 / n as f64);
        let b = h.point(h.t_max() * (k + 1) as f64 / n as f64);
        out.push((a - com, b - com, h.wire));
    }
    let end = h.point(h.t_max());
    let base = V3::new(0.0, 0.0, sw.head_z0() + sw.head_radius);
    let mid = (end + base) * 0.5;
    out.push((end - com, mid - com, h.wire));
    out.push((mid - com, base - com, h.wire));
    let top = V3::new(0.0, 0.0, sw.nose_z() - sw.head_radius);
    out.push((base - com, top - com, sw.head_radius));
    out
}

/// The contact shape sim-core sees: spheres along each capsule, no farther apart than their
/// radius (the surface dips under 0.14 radius between them). Spheres because sim-core 0.9's
/// capsule-box contact pulls a capsule into the box, and capsules fall through a cylinder's flat
/// end (`steer contact`; docs/CORTENFORGE.md); sphere contacts work on every shape.
fn spheres(caps: &[(V3, V3, f64)]) -> Vec<(V3, f64)> {
    let mut out: Vec<(V3, f64)> = Vec::new();
    for (a, c, r) in caps {
        let n = ((c - a).norm() / r).ceil().max(1.0) as usize;
        for i in 0..=n {
            let p = a + (c - a) * (i as f64 / n as f64);
            if out.last().is_none_or(|(q, rq)| (p - q).norm() > 1e-3 * r || rq != r) {
                out.push((p, *r));
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------------------------
// The field.

/// The coils' field per amp on a grid over the glycerol, with its gradient, for each pair.
pub struct Map {
    lo: V3,
    h: f64,
    n: usize,
    /// Per node: for each pair, the field per amp and its Jacobian (J_ij = dB_i / dx_j).
    b: Vec<[(V3, Matrix3<f64>); 3]>,
    /// Each pair's centre field per amp, along its axis.
    k: [f64; 3],
}

impl Map {
    pub fn new(h: f64) -> Map {
        let pairs = coils::pairs();
        let half = LIQ_R.max(LIQ_H) + h;
        let n = (2.0 * half / h).round() as usize + 1;
        let lo = V3::new(-half, -half, -half);
        let dl = 0.25e-3;
        let nodes: Vec<V3> = (0..n * n * n)
            .map(|idx| {
                let (i, j, l) = (idx / (n * n), (idx / n) % n, idx % n);
                lo + V3::new(i as f64, j as f64, l as f64) * h
            })
            .collect();
        let per_pair: Vec<Vec<(V3, Matrix3<f64>)>> = pairs
            .iter()
            .map(|p| {
                nodes
                    .par_chunks(256)
                    .flat_map_iter(|chunk| {
                        let mut pts = Vec::with_capacity(chunk.len() * 7);
                        for q in chunk {
                            pts.push(*q);
                            for ax in 0..3 {
                                let mut e = V3::zeros();
                                e[ax] = dl;
                                pts.push(q + e);
                                pts.push(q - e);
                            }
                        }
                        let b = p.field(&pts);
                        (0..chunk.len())
                            .map(|c| {
                                let s = 7 * c;
                                let mut j = Matrix3::zeros();
                                for ax in 0..3 {
                                    j.set_column(ax, &((b[s + 1 + 2 * ax] - b[s + 2 + 2 * ax]) / (2.0 * dl)));
                                }
                                (b[s], j)
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect()
            })
            .collect();
        let b = (0..nodes.len())
            .map(|i| [per_pair[0][i], per_pair[1][i], per_pair[2][i]])
            .collect();
        let k = [0, 1, 2].map(|i| pairs[i].field(&[V3::zeros()])[0].dot(&pairs[i].axis));
        Map { lo, h, n, b, k }
    }

    /// Trilinear interpolation of each pair's field and Jacobian per amp at x.
    fn at(&self, x: &V3) -> [(V3, Matrix3<f64>); 3] {
        let n = self.n;
        let mut i0 = [0usize; 3];
        let mut t = [0.0; 3];
        for a in 0..3 {
            let u = ((x[a] - self.lo[a]) / self.h).clamp(0.0, (n - 1) as f64 - 1e-9);
            i0[a] = u.floor() as usize;
            t[a] = u - i0[a] as f64;
        }
        let mut out = [(V3::zeros(), Matrix3::zeros()); 3];
        for c in 0..8 {
            let (di, dj, dk) = (c >> 2 & 1, c >> 1 & 1, c & 1);
            let w = (if di == 1 { t[0] } else { 1.0 - t[0] })
                * (if dj == 1 { t[1] } else { 1.0 - t[1] })
                * (if dk == 1 { t[2] } else { 1.0 - t[2] });
            let node = &self.b[((i0[0] + di) * n + i0[1] + dj) * n + i0[2] + dk];
            for p in 0..3 {
                out[p].0 += node[p].0 * w;
                out[p].1 += node[p].1 * w;
            }
        }
        out
    }
}

/// What makes the field: a uniform one, or the coils with each pair's current set so the centre
/// field is the one asked for.
pub enum Field<'a> {
    Uniform,
    Coils(&'a Map),
}

impl Field<'_> {
    /// The field and its Jacobian at x when the centre field is b.
    fn at(&self, x: &V3, b: &V3) -> (V3, Matrix3<f64>) {
        match self {
            Field::Uniform => (*b, Matrix3::zeros()),
            Field::Coils(map) => {
                let per = map.at(x);
                let mut bf = V3::zeros();
                let mut j = Matrix3::zeros();
                for p in 0..3 {
                    let amps = b[p] / map.k[p];
                    bf += per[p].0 * amps;
                    j += per[p].1 * amps;
                }
                (bf, j)
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// The simulation.

/// The swimmer in sim-core. `env` is MJCF for the world's fixed geoms (contype 0, conaffinity
/// 1); the swimmer's geoms collide with them only.
pub struct Sim<'a> {
    model: Model,
    data: Data,
    b: &'a Body,
    field: &'a Field<'a>,
    /// Gravity (m/s^2; 0 for a weightless swimmer, buoyancy off too).
    g: f64,
    dt: f64,
    k: Matrix6<f64>,
    mm: Matrix6<f64>,
    /// The field's direction now, the way the swimmer should go (the field turns about it,
    /// right-handed), its size (T) and rate (Hz).
    pub bhat: V3,
    pub n: V3,
    pub b0: f64,
    pub f: f64,
    pub t: f64,
    /// Largest relative miss between a contact-free step's end velocity and the prediction.
    pub miss: f64,
    /// Steps with a contact, of all steps.
    pub contact_steps: usize,
    pub steps: usize,
    /// Plain explicit drag instead of the linearly implicit one (needs steps of ~20 us).
    pub explicit: bool,
}

fn quat_from_axis(axis: &V3, x: &V3) -> UnitQuaternion<f64> {
    let z = axis.normalize();
    let xb = (x - z * z.dot(x)).normalize();
    let rot = Rotation3::from_matrix_unchecked(Matrix3::from_columns(&[xb, z.cross(&xb), z]));
    UnitQuaternion::from_rotation_matrix(&rot)
}

impl<'a> Sim<'a> {
    /// The swimmer with its centre of mass at `pos`, axis along `axis`, at rest, its moment
    /// lined up with the field (which starts along `x`, made square to `axis`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        b: &'a Body,
        field: &'a Field<'a>,
        dt: f64,
        env: &str,
        pos: V3,
        axis: V3,
        x: V3,
        friction: f64,
        gravity: bool,
        solref: f64,
    ) -> Sim<'a> {
        let g = if gravity { G } else { 0.0 };
        let q = quat_from_axis(&axis, &x);
        let i = &b.inertia;
        let mut geoms = String::new();
        for (p, r) in &b.spheres {
            geoms += &format!(
                "      <geom type=\"sphere\" pos=\"{} {} {}\" size=\"{r}\" contype=\"1\" conaffinity=\"0\" friction=\"{friction} 0.005 0.0001\" solref=\"{solref} 1\"/>\n",
                p.x, p.y, p.z
            );
        }
        let xml = format!(
            r#"<mujoco model="steer">
  <option timestep="{dt}" gravity="0 0 {}" viscosity="0" density="0"/>
  <worldbody>
{env}    <body name="swimmer" pos="{} {} {}" quat="{} {} {} {}">
      <freejoint/>
      <inertial pos="0 0 0" mass="{}" fullinertia="{} {} {} {} {} {}"/>
{geoms}    </body>
  </worldbody>
</mujoco>"#,
            -g,
            pos.x,
            pos.y,
            pos.z,
            q.w,
            q.i,
            q.j,
            q.k,
            b.mass,
            i[(0, 0)],
            i[(1, 1)],
            i[(2, 2)],
            i[(0, 1)],
            i[(0, 2)],
            i[(1, 2)],
        );
        let model = load_model(&xml).expect("steer model loads");
        let mut data = model.make_data();
        data.forward(&model).expect("forward");
        let mut mm = Matrix6::zeros();
        mm.fixed_view_mut::<3, 3>(0, 0).copy_from(&(Matrix3::identity() * b.mass));
        mm.fixed_view_mut::<3, 3>(3, 3).copy_from(&b.inertia);
        let k = (mm + b.r * dt).try_inverse().expect("M + dt R invertible");
        let bhat = (q * V3::x()).normalize();
        Sim {
            model,
            data,
            b,
            field,
            g,
            dt,
            k,
            mm,
            bhat,
            n: axis.normalize(),
            b0: 0.0,
            f: 0.0,
            t: 0.0,
            miss: 0.0,
            contact_steps: 0,
            steps: 0,
            explicit: false,
        }
    }

    /// Sim-core's own mass matrix against ours (largest relative difference).
    pub fn mass_matrix_check(&self) -> f64 {
        let qm = &self.data.qM;
        let mut worst = 0.0_f64;
        let scale = self.mm.abs().max();
        for r in 0..6 {
            for c in 0..6 {
                let ours = self.mm[(r, c)];
                worst = worst.max((qm[(r, c)] - ours).abs() / scale.max(ours.abs()).max(1e-300));
            }
        }
        worst
    }

    pub fn pos(&self) -> V3 {
        V3::new(self.data.qpos[0], self.data.qpos[1], self.data.qpos[2])
    }
    pub fn quat(&self) -> UnitQuaternion<f64> {
        let q = &self.data.qpos;
        UnitQuaternion::from_quaternion(Quaternion::new(q[3], q[4], q[5], q[6]))
    }
    pub fn axis(&self) -> V3 {
        self.quat() * V3::z()
    }
    /// The world velocity of the centre of mass and the body-frame spin.
    pub fn vel(&self) -> (V3, V3) {
        let v = &self.data.qvel;
        (V3::new(v[0], v[1], v[2]), V3::new(v[3], v[4], v[5]))
    }
    pub fn ncon(&self) -> usize {
        self.data.ncon
    }
    pub fn contacts(&self) -> &[cortenforge::sim::core::Contact] {
        &self.data.contacts[..self.data.ncon]
    }

    pub fn step(&mut self) {
        let b = self.b;
        let q = self.quat();
        let qi = q.inverse();
        let x = self.pos();
        let (v_w, w_b) = self.vel();
        let r_mag = q * b.mag;
        let (bf, j) = self.field.at(&(x + r_mag), &(self.bhat * self.b0));
        let m_w = q * V3::x() * b.m;
        let f_mag = j.transpose() * m_w;
        let t_mag = m_w.cross(&bf) + r_mag.cross(&f_mag);
        let f_buoy = V3::new(0.0, 0.0, GLYCEROL_RHO * b.volume * self.g);
        let t_buoy = (q * b.cob).cross(&f_buoy);
        let grav = V3::new(0.0, 0.0, -b.mass * self.g);
        let fk = qi * (f_mag + f_buoy + grav);
        let tk = qi * (t_mag + t_buoy);
        let v_b = qi * v_w;
        let s = Vector6::new(v_b.x, v_b.y, v_b.z, w_b.x, w_b.y, w_b.z);
        let kn = Vector6::new(fk.x, fk.y, fk.z, tk.x, tk.y, tk.z);
        let y = self.k * (self.mm * s + kn * self.dt);
        // Explicit (plain) drag on the velocity now, for checking at very short steps.
        let d = if self.explicit { -(b.r * s) } else { -(b.r * y) };
        let fd = q * V3::new(d[0], d[1], d[2]);
        let td = q * V3::new(d[3], d[4], d[5]);
        let t = t_mag + t_buoy + td;
        let f = f_mag + f_buoy + fd;
        self.data.xfrc_applied[1] = Vector6::new(t.x, t.y, t.z, f.x, f.y, f.z);
        self.data.step(&self.model).expect("step");
        self.steps += 1;
        if self.data.ncon > 0 {
            self.contact_steps += 1;
        } else {
            let (v2, w2) = self.vel();
            let v2b = qi * v2;
            let got = Vector6::new(v2b.x, v2b.y, v2b.z, w2.x, w2.y, w2.z);
            let scale = y.norm().max(1e-300);
            self.miss = self.miss.max((got - y).norm() / scale);
        }
        // Turn the field about n, kept square to n (so a change of n never jumps the field).
        let a = 2.0 * PI * self.f * self.dt;
        let n = self.n;
        let bh = self.bhat * a.cos() + n.cross(&self.bhat) * a.sin();
        let bh = bh - n * n.dot(&bh);
        if bh.norm() > 1e-9 {
            self.bhat = bh.normalize();
        } else {
            // n turned onto the field's direction: pick any direction square to n.
            self.bhat = coils::plane(&n).0;
        }
        self.t += self.dt;
    }

    /// Run for `secs`.
    pub fn run(&mut self, secs: f64) {
        let n = (secs / self.dt).round() as usize;
        for _ in 0..n {
            self.step();
        }
    }
}


// ---------------------------------------------------------------------------------------------

fn deg(a: &V3, b: &V3) -> f64 {
    a.normalize().dot(&b.normalize()).clamp(-1.0, 1.0).acos().to_degrees()
}

/// Means over a stretch of a run: the centre of mass's velocity (world), the spin about the
/// axis (Hz), and the axis's angle from a direction (mean and largest, deg).
struct Stats {
    v: V3,
    spin: f64,
    off: f64,
    off_max: f64,
}

fn measure(s: &mut Sim, span: f64, dir: &V3) -> Stats {
    let n = (span / s.dt).round().max(1.0) as usize;
    let x0 = s.pos();
    let (mut spin, mut off, mut off_max) = (0.0, 0.0, 0.0_f64);
    for _ in 0..n {
        s.step();
        let a = s.axis();
        let o = deg(&a, dir);
        spin += s.vel().1.z;
        off += o;
        off_max = off_max.max(o);
    }
    let k = n as f64;
    Stats {
        v: (s.pos() - x0) / (k * s.dt),
        spin: spin / k / (2.0 * PI),
        off: off / k,
        off_max,
    }
}

/// Free fluid: no walls.
fn free<'a>(b: &'a Body, field: &'a Field<'a>, dt: f64, axis: V3, gravity: bool) -> Sim<'a> {
    let x = if axis.x.abs() < 0.9 { V3::x() } else { V3::y() };
    Sim::new(b, field, dt, "", V3::zeros(), axis, x, 0.2, gravity, 0.02)
}

/// Whole field turns lasting about `secs`.
fn turns_span(f: f64, secs: f64) -> f64 {
    (secs * f).round().max(1.0) / f
}


// ---------------------------------------------------------------------------------------------
// The course.

/// The printed insert: a post on the jar's axis and two fins from it to the wall (at +x and
/// -x), each with a window, the second higher: a spiral staircase. Metres.
const POST_R: f64 = 4e-3;
const FIN_T: f64 = 1.6e-3;
/// Windows: radial span, height, and the centres' heights (fin at +x, fin at -x).
const WIN_R: (f64, f64) = (9e-3, 23.5e-3);
const WIN_H: f64 = 18e-3;
const WIN_Z: [f64; 2] = [-10e-3, 10e-3];
/// The path's radius (the swimmer's centre of mass laps the post at this radius).
const LAP_R: f64 = 14e-3;
/// The path is level this far (deg) each side of a window.
const PATH_LEVEL: f64 = 25.0;
/// Boxes in the polygon that stands for the jar's glass wall.
const WALL_N: usize = 48;
/// Resin on glass and resin, wetted with glycerol (a guess; contact is rare on a good run).
const FRICTION: f64 = 0.2;
/// Contact time constant (s): sim-core's default. With sphere contacts the swimmer presses under
/// 0.04 mm into a wall at 1 ms steps (`steer wall`).
const SOLREF: f64 = 0.02;
/// The camera(s): frames per second, position and axis noise (1 sigma, m and rad).
const FPS: f64 = 30.0;
const CAM_POS: f64 = 0.05e-3;
const CAM_AXIS: f64 = 0.5 * PI / 180.0;
/// The controller's field axis stays within this of its aim (rad).
const MAX_LEAD: f64 = 30.0 * PI / 180.0;
/// The controller keeps the field's rate under this fraction of step-out: turning and the
/// lead of the field's axis add to the twist, and a swimmer that steps out tumbles.
const CTRL_KAPPA: f64 = 0.6;
/// Standing up off the floor: a still field turned back over the swimmer at this rate (Hz)
/// through this angle (rad), levering the nose up about the tail.
const STAND_RATE: f64 = 0.05;
/// A run stops if it makes no progress along the path for this long (s).
const STALL: f64 = 15.0;
const STAND_ANGLE: f64 = 60.0 * PI / 180.0;

fn box_sdf(p: &V3, c: &V3, hs: &V3) -> f64 {
    let q = (p - c).abs() - hs;
    let out = V3::new(q.x.max(0.0), q.y.max(0.0), q.z.max(0.0)).norm();
    out + q.x.max(q.y).max(q.z).min(0.0)
}

#[derive(Clone)]
struct Course {
    boxes: Vec<(V3, V3)>,
    path: Vec<V3>,
    cum: Vec<f64>,
}

/// What a point is closest to.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Near {
    Insert,
    Wall,
    Floor,
    Surface,
}

impl Course {
    fn new(z_start: f64) -> Course {
        let mut boxes = Vec::new();
        let hh = LIQ_H + 1e-3;
        let (x_in, x_out) = (POST_R - 0.5e-3, LIQ_R + 1e-3);
        for (side, zw) in [(1.0, WIN_Z[0]), (-1.0, WIN_Z[1])] {
            let (z_lo, z_hi) = (zw - 0.5 * WIN_H, zw + 0.5 * WIN_H);
            let mut add = |x0: f64, x1: f64, z0: f64, z1: f64| {
                boxes.push((
                    V3::new(side * 0.5 * (x0 + x1), 0.0, 0.5 * (z0 + z1)),
                    V3::new(0.5 * (x1 - x0), 0.5 * FIN_T, 0.5 * (z1 - z0)),
                ));
            };
            add(x_in, x_out, -hh, z_lo);
            add(x_in, x_out, z_hi, hh);
            add(x_in, WIN_R.0, z_lo, z_hi);
            add(WIN_R.1, x_out, z_lo, z_hi);
        }
        let mut c = Course {
            boxes,
            path: Vec::new(),
            cum: Vec::new(),
        };
        c.set_path(-90.0, z_start);
        c
    }

    /// The path from (phi0 deg, z0): anticlockwise round the post, level through the first
    /// window (0 deg) and the second (180 deg), climbing between them, to above the start
    /// (270 deg).
    fn set_path(&mut self, phi0: f64, z0: f64) {
        let knots = [
            (phi0, z0),
            (-PATH_LEVEL, WIN_Z[0]),
            (PATH_LEVEL, WIN_Z[0]),
            (180.0 - PATH_LEVEL, WIN_Z[1]),
            (180.0 + PATH_LEVEL, WIN_Z[1]),
            (270.0, WIN_Z[1]),
        ];
        let mut path = Vec::new();
        for w in knots.windows(2) {
            let ((p0, z0), (p1, z1)) = (w[0], w[1]);
            let n = ((p1 - p0).round() as usize).max(1);
            for i in 0..n {
                let t = i as f64 / n as f64;
                let phi = (p0 + t * (p1 - p0)).to_radians();
                path.push(V3::new(LAP_R * phi.cos(), LAP_R * phi.sin(), z0 + t * (z1 - z0)));
            }
        }
        let (pl, zl) = knots[knots.len() - 1];
        path.push(V3::new(LAP_R * pl.to_radians().cos(), LAP_R * pl.to_radians().sin(), zl));
        let mut cum = vec![0.0];
        for w in path.windows(2) {
            cum.push(cum.last().unwrap() + (w[1] - w[0]).norm());
        }
        self.path = path;
        self.cum = cum;
    }

    fn length(&self) -> f64 {
        *self.cum.last().unwrap()
    }

    /// Distance from a point to the nearest solid (the insert, the glass, the floor, the
    /// glycerol's surface), and which.
    fn clearance(&self, p: &V3) -> (f64, Near) {
        let r = (p.x * p.x + p.y * p.y).sqrt();
        let mut best = (r - POST_R, Near::Insert);
        for (c, hs) in &self.boxes {
            let d = box_sdf(p, c, hs);
            if d < best.0 {
                best = (d, Near::Insert);
            }
        }
        for (d, w) in [
            (LIQ_R - r, Near::Wall),
            (p.z + LIQ_H, Near::Floor),
            (LIQ_H - p.z, Near::Surface),
        ] {
            if d < best.0 {
                best = (d, w);
            }
        }
        best
    }

    /// MJCF for the jar and the insert; returns it and the number of geoms (they come first).
    fn mjcf(&self) -> (String, usize) {
        let fr = format!(
            "contype=\"0\" conaffinity=\"1\" friction=\"{FRICTION} 0.005 0.0001\" solref=\"{SOLREF} 1\""
        );
        let mut s = String::new();
        let mut n = 0;
        s += &format!(
            "    <geom type=\"plane\" pos=\"0 0 {}\" size=\"0.04 0.04 0.001\" {fr}/>\n",
            -LIQ_H
        );
        s += &format!(
            "    <geom type=\"box\" pos=\"0 0 {}\" size=\"0.04 0.04 0.001\" {fr}/>\n",
            LIQ_H + 1e-3
        );
        n += 2;
        let w = 1e-3;
        let half_t = LIQ_R * (PI / WALL_N as f64).tan() * 1.02;
        for i in 0..WALL_N {
            let a = 2.0 * PI * i as f64 / WALL_N as f64;
            let rc = LIQ_R + w;
            s += &format!(
                "    <geom type=\"box\" pos=\"{} {} 0\" quat=\"{} 0 0 {}\" size=\"{w} {half_t} {}\" {fr}/>\n",
                rc * a.cos(),
                rc * a.sin(),
                (0.5 * a).cos(),
                (0.5 * a).sin(),
                LIQ_H + 1e-3
            );
            n += 1;
        }
        s += &format!(
            "    <geom type=\"cylinder\" pos=\"0 0 0\" size=\"{POST_R} {}\" {fr}/>\n",
            LIQ_H + 1e-3
        );
        n += 1;
        for (c, hs) in &self.boxes {
            s += &format!(
                "    <geom type=\"box\" pos=\"{} {} {}\" size=\"{} {} {}\" {fr}/>\n",
                c.x, c.y, c.z, hs.x, hs.y, hs.z
            );
            n += 1;
        }
        (s, n)
    }
}

/// The coil set's numbers for the drive: centre field per amp, resistance, copper mass.
struct Drive {
    k: [f64; 3],
    r: [f64; 3],
    cu: [f64; 3],
}

impl Drive {
    fn new() -> Drive {
        let p = coils::pairs();
        Drive {
            k: [0, 1, 2].map(|i| p[i].field(&[V3::zeros()])[0].dot(&p[i].axis)),
            r: [0, 1, 2].map(|i| p[i].coil.resistance()),
            cu: [0, 1, 2].map(|i| p[i].coil.copper_mass()),
        }
    }
    /// Power in each pair for centre field b (W).
    fn power(&self, b: &V3) -> [f64; 3] {
        [0, 1, 2].map(|i| self.r[i] * (b[i] / self.k[i]).powi(2))
    }
}

/// A small deterministic normal sampler (xorshift64* and Box-Muller).
struct Rng(u64);

impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        ((self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        let (u, v) = (self.uniform(), self.uniform());
        (-2.0 * u.ln()).sqrt() * (2.0 * PI * v).cos()
    }
    fn vec(&mut self, s: f64) -> V3 {
        V3::new(self.normal(), self.normal(), self.normal()) * s
    }
}

/// The controller, once a camera frame: pure pursuit along the path, then (if `comp`) the
/// swimmer's own sinking taken out of the aim (from the model's drag and its seen axis) and the
/// axis's offset from the field's (droop, swing, lag) taken out by a PI loop on the seen axis.
/// Without `comp` it points the field's axis at the target.
#[derive(Clone, Copy)]
struct Gains {
    /// (m/s) the speed asked for along the path; (m) how far ahead the target is.
    v: f64,
    look: f64,
    /// The axis loop: proportional, and integral (1/s).
    kp: f64,
    ki: f64,
    /// Vertical integral on the height error (1/s).
    kz: f64,
    comp: bool,
    /// Stand up off the floor with a still field first.
    stand: bool,
}

struct Ctrl {
    g: Gains,
    idx: usize,
    corr: V3,
    zbias: f64,
    mu_par: f64,
    mu_perp: f64,
}

impl Ctrl {
    fn new(g: Gains, b: &Body) -> Ctrl {
        let a_inv = b
            .r
            .fixed_view::<3, 3>(0, 0)
            .into_owned()
            .try_inverse()
            .expect("drag invertible");
        Ctrl {
            g,
            idx: 0,
            corr: V3::zeros(),
            zbias: 0.0,
            mu_par: a_inv[(2, 2)],
            mu_perp: 0.5 * (a_inv[(0, 0)] + a_inv[(1, 1)]),
        }
    }

    /// The closest path index near the last one (progress only goes forward).
    fn progress(&mut self, c: &Course, x: &V3) -> usize {
        let hi = (self.idx + 40).min(c.path.len() - 1);
        let mut best = (f64::INFINITY, self.idx);
        for i in self.idx..=hi {
            let d = (c.path[i] - x).norm();
            if d < best.0 {
                best = (d, i);
            }
        }
        self.idx = best.1;
        best.1
    }

    /// (n, f) from the seen centre of mass and axis.
    fn update(&mut self, c: &Course, b: &Body, x: &V3, a: &V3, f_max: f64, dt: f64) -> (V3, f64) {
        let i = self.progress(c, x);
        let s_look = c.cum[i] + self.g.look;
        let j = c.cum.partition_point(|&s| s < s_look).min(c.path.len() - 1);
        let d = (c.path[j] - x).normalize();
        if !self.g.comp {
            return (d, (self.g.v / b.v_hz_free).min(f_max));
        }
        self.zbias += self.g.kz * dt * (c.path[i].z - x.z);
        self.zbias = self.zbias.clamp(-3e-3, 3e-3);
        let w = b.weight();
        let az = a.z;
        let sink = -(a * (self.mu_par * az) + (V3::z() - a * az) * self.mu_perp) * w;
        let v_cmd = d * self.g.v - sink + V3::z() * self.zbias;
        let f = (v_cmd.norm() / b.v_hz_free).min(f_max);
        let a_des = v_cmd.normalize();
        // The swimmer moves along the field's axis nose first or tail first (its magnet turned
        // either way round): steer the way it actually goes.
        let a = if a.dot(&a_des) < 0.0 { -a } else { *a };
        let e = a_des - a;
        let e = e - a_des * a_des.dot(&e);
        self.corr += e * (self.g.ki * dt);
        self.corr -= a_des * a_des.dot(&self.corr);
        if self.corr.norm() > 1.0 {
            self.corr = self.corr.normalize();
        }
        // The field's axis: the aim plus the loop's correction, at most MAX_LEAD from the aim.
        let n = (a_des + e * self.g.kp + self.corr).normalize();
        let lead = a_des.dot(&n).clamp(-1.0, 1.0).acos();
        if lead > MAX_LEAD {
            let side = (n - a_des * a_des.dot(&n)).normalize();
            return ((a_des * MAX_LEAD.cos() + side * MAX_LEAD.sin()).normalize(), f);
        }
        (n, f)
    }
}

/// One run of the course.
struct Outcome {
    finished: bool,
    time: f64,
    err_rms: f64,
    err_max: f64,
    /// Closest the body came to the insert, the glass, and the floor after lifting off (m).
    near_insert: f64,
    near_wall: f64,
    /// Steps touching after lift-off, and the solids touched.
    touch_secs: f64,
    touched: Vec<&'static str>,
    watts: f64,
    heat: f64,
    progress: f64,
    wall_secs: f64,
    trace: Vec<(f64, V3, V3, f64)>,
    /// The axis's elevation after standing up (rad), and when the body left the floor (s).
    stand_elev: f64,
    t_lift: f64,
}

#[allow(clippy::too_many_arguments)]
fn course_run(
    b: &Body,
    field: &Field,
    c: &Course,
    g: Gains,
    b0: f64,
    dt: f64,
    seed: u64,
    t_max: f64,
) -> Outcome {
    let t_wall = Instant::now();
    let (env, n_env) = c.mjcf();
    // The path is for the swimmer's middle (the centre of mass is near the nose).
    let start = c.path[0] - V3::x() * b.mid;
    let mut s = Sim::new(b, field, dt, &env, start + V3::new(0.0, 0.0, 0.2e-3), V3::x(), V3::z(), FRICTION, true, SOLREF);
    // Settle onto the floor with the field off.
    s.run(1.0);
    let t0 = s.t;
    let drive = Drive::new();
    let mut energy = [0.0; 3];
    let advance = |s: &mut Sim, secs: f64, energy: &mut [f64; 3]| {
        for _ in 0..(secs / dt).round() as usize {
            s.step();
            let p = drive.power(&(s.bhat * s.b0));
            for k in 0..3 {
                energy[k] += p[k] * dt;
            }
        }
    };
    // Stand up: a still field pointing up rolls the magnet to the top, then turns slowly back
    // over the swimmer about the horizontal line across it, and the magnet levers the nose up
    // about the tail. (Without it the turning field drags the swimmer along the floor and the
    // insert before it gets clear: the `stand: false` run.)
    let a0 = s.axis();
    let t_hat = V3::new(a0.x, a0.y, 0.0).normalize();
    s.b0 = b0;
    s.n = t_hat;
    if g.stand {
        s.n = t_hat.cross(&V3::z()).normalize();
        s.bhat = V3::z();
        s.f = 0.0;
        advance(&mut s, 0.5, &mut energy);
        s.f = STAND_RATE;
        advance(&mut s, STAND_ANGLE / (2.0 * PI * STAND_RATE), &mut energy);
        s.f = 0.0;
        advance(&mut s, 1.0, &mut energy);
    }
    let stand_elev = s.axis().z.asin();
    // The path starts where the stand-up left the swimmer's middle.
    let mid0 = s.pos() + s.axis() * b.mid;
    let mut course = c.clone();
    course.set_path(mid0.y.atan2(mid0.x).to_degrees(), mid0.z);
    let c = &course;
    let mut ctrl = Ctrl::new(g, b);
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15 ^ seed);
    let f_max = CTRL_KAPPA * b.fso_t_free * b0;
    let frame = 1.0 / FPS;
    let mut next_frame = s.t;
    let mut pending: Option<(V3, f64)> = None;
    let mut lifted = false;
    let mut t_lift = f64::NAN;
    let (mut err2, mut err_max, mut nerr) = (0.0, 0.0_f64, 0usize);
    let (mut near_insert, mut near_wall) = (f64::INFINITY, f64::INFINITY);
    let mut touch_steps = 0usize;
    let mut touched: Vec<&'static str> = Vec::new();
    let mut finished = false;
    let mut trace = Vec::new();
    let end = *c.path.last().unwrap();
    let (mut best_idx, mut t_best) = (0usize, s.t);
    while s.t - t0 < t_max {
        if s.t >= next_frame {
            next_frame += frame;
            if let Some((n, f)) = pending {
                s.n = n;
                s.f = f;
            }
            let com = s.pos();
            let xt = com + s.axis() * b.mid;
            let x = xt + rng.vec(CAM_POS);
            let a = (s.axis() + rng.vec(CAM_AXIS)).normalize();
            pending = Some(ctrl.update(c, b, &x, &a, f_max, frame));
            // Path error and clearances, from the true state.
            let e = c.path.iter().map(|p| (p - xt).norm()).fold(f64::INFINITY, f64::min);
            err2 += e * e;
            err_max = err_max.max(e);
            nerr += 1;
            let q = s.quat();
            let mut floor = f64::INFINITY;
            for (pa, pc, rad) in &b.capsules {
                for t in [0.0, 0.5, 1.0] {
                    let p = com + q * (pa + (pc - pa) * t);
                    floor = floor.min(p.z + LIQ_H - rad);
                    let (d, w) = c.clearance(&p);
                    let d = d - rad;
                    match w {
                        Near::Insert => near_insert = near_insert.min(d),
                        Near::Wall => near_wall = near_wall.min(d),
                        _ => {}
                    }
                }
            }
            if !lifted && floor > 0.5e-3 {
                lifted = true;
                t_lift = s.t - t0;
            }
            trace.push((s.t - t0, xt, s.axis(), s.f));
            if ctrl.idx > best_idx {
                best_idx = ctrl.idx;
                t_best = s.t;
            } else if s.t - t_best > STALL {
                break;
            }
            if ctrl.idx + 1 >= c.path.len() - 1 || (xt - end).norm() < 2e-3 {
                finished = true;
                break;
            }
        }
        s.step();
        let p = drive.power(&(s.bhat * s.b0));
        for k in 0..3 {
            energy[k] += p[k] * dt;
        }
        if lifted && s.ncon() > 0 {
            touch_steps += 1;
            for con in s.contacts() {
                let gid = con.geom1.min(con.geom2);
                let what = if gid < n_env {
                    match gid {
                        0 => "floor",
                        1 => "surface",
                        g if g < 2 + WALL_N => "glass",
                        g if g == 2 + WALL_N => "post",
                        _ => "fin",
                    }
                } else {
                    "itself"
                };
                if !touched.contains(&what) {
                    touched.push(what);
                }
            }
        }
    }
    let time = s.t - t0;
    let joules: f64 = energy.iter().sum();
    let heat = (0..3)
        .map(|k| energy[k] / (drive.cu[k] * 385.0))
        .fold(0.0_f64, f64::max);
    Outcome {
        finished,
        time,
        err_rms: (err2 / nerr.max(1) as f64).sqrt(),
        err_max,
        near_insert,
        near_wall,
        touch_secs: touch_steps as f64 * dt,
        touched,
        watts: joules / time,
        heat,
        progress: c.cum[ctrl.idx] / c.length(),
        wall_secs: t_wall.elapsed().as_secs_f64(),
        trace,
        stand_elev,
        t_lift,
    }
}

fn show(label: &str, o: &Outcome) {
    let touched = if o.touched.is_empty() {
        "-".to_string()
    } else {
        o.touched.join("+")
    };
    println!(
        "   {label:<34} {:>6.1} {:>6.1} {:>8} {:>7.1} {:>11} {:>8.2} {:>8.2} {:>10.2} {:>7} {:>8.1} {:>8.2} {:>7.0}",
        o.stand_elev.to_degrees(),
        o.t_lift,
        if o.finished {
            "yes".to_string()
        } else {
            format!("{:.0}%", 100.0 * o.progress)
        },
        o.time,
        format!("{:.2} / {:.2}", o.err_rms * 1e3, o.err_max * 1e3),
        o.near_insert * 1e3,
        o.near_wall * 1e3,
        o.touch_secs,
        touched,
        o.watts,
        o.heat,
        o.wall_secs
    );
}

fn header() {
    println!(
        "   {:<34} {:>6} {:>6} {:>8} {:>7} {:>11} {:>8} {:>8} {:>10} {:>7} {:>8} {:>8} {:>7}",
        "run", "stood", "off s", "finished", "time s", "path mm", "insert", "glass", "touch s", "touched", "mean W", "worst C", "wall s"
    );
}

pub fn run() {
    let t_all = Instant::now();
    let sw = shape::chosen();
    let b = Body::new(sw);
    let h = sw.helix;
    println!(
        "steer: the swimmer `shape` chose (R {:.2} mm, pitch {:.2} mm, wire {:.2} mm, {} turns, {:.1} mm, {:.1} mg) in sim-core",
        h.radius * 1e3,
        h.pitch * 1e3,
        2.0 * h.wire * 1e3,
        h.turns,
        sw.length() * 1e3,
        b.mass * 1e6
    );
    println!(
        "   model: {} points at 16 round, extrapolated from 12 and 16 ({:.0} s)",
        b.model_points, b.model_secs
    );
    println!(
        "   axis held: {:.3} mm/s per Hz, step-out {:.3} Hz per mT; axis free: {:.3} mm/s per Hz, {:.3} Hz per mT",
        b.v_hz * 1e3,
        b.fso_t * 1e-3,
        b.v_hz_free * 1e3,
        b.fso_t_free * 1e-3
    );
    println!(
        "   sinking upright {:.2} mm/s, level {:.2} mm/s (orientation held); net weight {:.1} mg in glycerol",
        b.sink_up * 1e3,
        b.sink_level * 1e3,
        b.weight() / G * 1e6
    );
    println!(
        "   inertia about the centre of mass (g mm^2): {:.4} {:.4} {:.4}; sampled mass {:.2} mg against {:.2} ({:+.2}%)",
        b.inertia[(0, 0)] * 1e9,
        b.inertia[(1, 1)] * 1e9,
        b.inertia[(2, 2)] * 1e9,
        b.mass_sampled * 1e6,
        b.mass * 1e6,
        100.0 * (b.mass_sampled / b.mass - 1.0)
    );
    let a_tt = b.r.fixed_view::<3, 3>(0, 0).into_owned();
    let tau = b.mass / a_tt.symmetric_eigenvalues().max();
    let tau_r = b.inertia[(2, 2)] / b.r[(5, 5)];
    println!(
        "   it settles to a new speed in {:.2} ms and a new spin in {:.3} ms (mass / drag)",
        tau * 1e3,
        tau_r * 1e3,
    );
    let (tail, nose) = b.ends();
    // Sideways drag: unit velocity across the axis; the torque about the centre of mass over the
    // force is how far behind it the drag acts.
    let lever = -b.r[(4, 0)] / b.r[(0, 0)];
    println!(
        "   from the centre of mass: tail {:.2} mm, nose {:+.2} mm, middle {:.2} mm, buoyancy's centre {:.2} mm, magnet {:+.2} mm; sideways drag acts {:.2} mm behind it",
        tail * 1e3,
        nose * 1e3,
        b.mid * 1e3,
        b.cob.z * 1e3,
        b.mag.z * 1e3,
        lever * 1e3
    );
    let uni = Field::Uniform;
    println!(
        "   sim-core's mass matrix against ours: {:.1e}",
        free(&b, &uni, 1e-3, V3::z(), true).mass_matrix_check()
    );

    // -----------------------------------------------------------------------------------------
    let dt = 1e-3;
    // `steer course`: the course alone (sections 1 and 2 take about 30 s).
    if std::env::args().nth(2).as_deref() != Some("course") {
    let b2 = 2e-3;
    let f_so = b.fso_t_free * b2;
    let f = KAPPA * f_so;
    println!(
        "1. The coupling, checked (uniform 2 mT field turning at {KAPPA} of the free axis's step-out, {f:.3} Hz)"
    );
    println!(
        "   swimming alone (weightless, axis up): the model {:.4} mm/s axis held, {:.4} free",
        b.v_hz * f * 1e3,
        b.v_hz_free * f * 1e3
    );
    println!(
        "   {:>10} {:>12} {:>10} {:>10} {:>10} {:>14} {:>12}",
        "step ms", "speed mm/s", "vs free", "spin Hz", "tilt deg", "wall s per s", "step miss"
    );
    for dt in [2e-5, 2e-4, 1e-3, 4e-3] {
        let t0 = Instant::now();
        let mut s = free(&b, &uni, dt, V3::z(), false);
        s.b0 = b2;
        s.f = f;
        s.run(turns_span(f, 1.0));
        let st = measure(&mut s, turns_span(f, 3.0), &V3::z());
        let secs = t0.elapsed().as_secs_f64() / s.t;
        println!(
            "   {:>10.3} {:>12.4} {:>+9.2}% {:>10.4} {:>10.3} {:>14.4} {:>12.1e}",
            dt * 1e3,
            st.v.z * 1e3,
            100.0 * (st.v.z / (b.v_hz_free * f) - 1.0),
            st.spin,
            st.off,
            secs,
            s.miss
        );
    }
    {
        let mut s = free(&b, &uni, dt, V3::x(), true);
        let mut t_down = f64::NAN;
        while s.t < 60.0 {
            s.step();
            if deg(&s.axis(), &-V3::z()) < 5.0 {
                t_down = s.t;
                break;
            }
        }
        s.run(5.0);
        let st = measure(&mut s, 5.0, &-V3::z());
        // Free to turn as it sinks (a chiral body spins as it falls): the full 6 x 6 mobility.
        let w = b.weight();
        let r_inv = b.r.try_inverse().expect("resistance invertible");
        let y = r_inv * Vector6::new(0.0, 0.0, -w, 0.0, 0.0, 0.0);
        println!(
            "   sinking alone (no field, released level): turns nose down (within 5 deg) in {t_down:.1} s, then sinks at {:.4} mm/s; the model {:.4} free to turn ({:.4} held), spinning {:.4} Hz (sim-core {:.4})",
            -st.v.z * 1e3,
            -y[2] * 1e3,
            b.sink_up * 1e3,
            y[5].abs() / (2.0 * PI),
            st.spin.abs()
        );
    }
    let mut s = free(&b, &uni, dt, V3::z(), true);
    s.b0 = b2;
    s.f = f;
    s.run(turns_span(f, 2.0));
    let st = measure(&mut s, turns_span(f, 4.0), &V3::z());
    println!(
        "   climbing (both): {:.4} mm/s; the model, swimming less sinking, {:.4} free ({:.4} held); axis {:.2} deg off vertical",
        st.v.z * 1e3,
        (b.v_hz_free * f - b.sink_up) * 1e3,
        (b.v_hz * f - b.sink_up) * 1e3,
        st.off
    );
    for k in [1.25, 1.5, 2.0] {
        let f_hi = k * f_so;
        let spin = crate::swim::spin(f_hi, f_so);
        let mut s = free(&b, &uni, dt, V3::z(), false);
        s.b0 = b2;
        s.f = f_hi;
        s.run(4.0);
        let st = measure(&mut s, turns_span(f_hi, 8.0), &V3::z());
        println!(
            "   past step-out, weightless ({k}x, {f_hi:.2} Hz): Adler's spin {spin:.3} Hz, sim-core {:.3} Hz; speed {:.3} mm/s against {:.3}; axis {:.1} deg off (worst {:.1})",
            st.spin,
            st.v.z * 1e3,
            b.v_hz_free * spin * 1e3,
            st.off,
            st.off_max
        );
    }

    // -----------------------------------------------------------------------------------------
    println!(
        "2. Swimming level and turning (free fluid; the field turning at {KAPPA} of step-out)"
    );
    println!(
        "   level along x, settled 30 s, then 10 s: where the axis points against the field's (nose down, swung sideways, in all)"
    );
    println!(
        "   {:>5} {:>8} {:>10} {:>10} {:>9} {:>12} {:>12}",
        "B mT", "f Hz", "down deg", "side deg", "all deg", "fwd mm/s", "sink mm/s"
    );
    let fields = [2e-3, 3e-3, 4e-3, 5e-3];
    for &b0 in &fields {
        let f = KAPPA * b.fso_t_free * b0;
        let mut s = free(&b, &uni, dt, V3::x(), true);
        s.b0 = b0;
        s.f = f;
        s.run(30.0);
        let n = (turns_span(f, 10.0) / dt).round() as usize;
        let x0 = s.pos();
        let (mut down, mut side, mut all) = (0.0, 0.0, 0.0);
        for _ in 0..n {
            s.step();
            let a = s.axis();
            down -= a.z.asin().to_degrees();
            side += a.y.atan2(a.x).to_degrees();
            all += deg(&a, &V3::x());
        }
        let k = n as f64;
        let v = (s.pos() - x0) / (k * dt);
        println!(
            "   {:>5.0} {:>8.2} {:>10.2} {:>10.2} {:>9.2} {:>12.3} {:>12.3}",
            b0 * 1e3,
            f,
            down / k,
            side / k,
            all / k,
            v.x * 1e3,
            -v.z * 1e3
        );
    }
    println!(
        "   (model with the axis held level: forward {:.3} mm/s per mT of field, sinking {:.3} mm/s)",
        b.v_hz_free * KAPPA * b.fso_t_free * 1e-3 * 1e3,
        b.sink_level * 1e3
    );
    println!(
        "   pitching up and turning: the field's axis swung at a steady rate from level along x (after 30 s level); kept up = within 60 deg of it throughout"
    );
    println!(
        "   {:>5} {:>12} {:>8} {:>9} {:>12} {:>12} {:>12} {:>10}",
        "B mT", "turn", "rad/s", "for rad", "behind deg", "worst deg", "radius mm", "kept up"
    );
    for &b0 in &fields {
        let f = KAPPA * b.fso_t_free * b0;
        let v = b.v_hz_free * f;
        for (name, kind, span) in [("pitch up", 1, 0.5 * PI), ("yaw", 0, PI)] {
            for om in [0.05, 0.1, 0.2, 0.3, 0.5] {
                let mut s = free(&b, &uni, dt, V3::x(), true);
                s.b0 = b0;
                s.f = f;
                s.run(30.0);
                let steps = (span / om / dt).round() as usize;
                let (mut behind, mut worst, mut cnt) = (0.0, 0.0_f64, 0usize);
                for i in 0..steps {
                    let a = om * (i as f64 + 1.0) * dt;
                    let (n, tdir) = match kind {
                        0 => (V3::new(a.cos(), a.sin(), 0.0), V3::new(-a.sin(), a.cos(), 0.0)),
                        _ => (V3::new(a.cos(), 0.0, a.sin()), V3::new(-a.sin(), 0.0, a.cos())),
                    };
                    s.n = n;
                    s.step();
                    let ax = s.axis();
                    worst = worst.max(deg(&ax, &n));
                    if i >= steps / 2 {
                        behind -= ax.dot(&tdir).clamp(-1.0, 1.0).asin().to_degrees();
                        cnt += 1;
                    }
                }
                let kept = worst < 60.0;
                println!(
                    "   {:>5.0} {:>12} {:>8.2} {:>9.2} {:>12.2} {:>12.2} {:>12.1} {:>10}",
                    b0 * 1e3,
                    name,
                    om,
                    span,
                    behind / cnt as f64,
                    worst,
                    v / om * 1e3,
                    if kept { "yes" } else { "no" }
                );
                if !kept {
                    break;
                }
            }
        }
    }
    println!("   ({:.0} s so far)", t_all.elapsed().as_secs_f64());
    }

    // -----------------------------------------------------------------------------------------
    let z_start = -LIQ_H + h.radius + h.wire;
    let c = Course::new(z_start);
    println!(
        "3. The course: a printed post ({:.0} mm across) with two fins ({:.1} mm thick) to the glass, each with a window {:.0} mm tall from {:.1} to {:.1} mm out",
        2.0 * POST_R * 1e3,
        FIN_T * 1e3,
        WIN_H * 1e3,
        WIN_R.0 * 1e3,
        WIN_R.1 * 1e3
    );
    println!(
        "   windows centred {:.0} mm (fin at +x) and {:+.0} mm (fin at -x) from mid-height; start lying on the floor; path {:.1} mm round the post at {:.0} mm",
        WIN_Z[0] * 1e3,
        WIN_Z[1] * 1e3,
        c.length() * 1e3,
        LAP_R * 1e3
    );
    println!(
        "   camera: {FPS:.0} frames/s, position noise {:.2} mm, axis noise {:.1} deg, one frame late; friction {FRICTION}",
        CAM_POS * 1e3,
        CAM_AXIS.to_degrees()
    );
    let climb = ((WIN_Z[1] - WIN_Z[0]) / (LAP_R * (180.0 - 2.0 * PATH_LEVEL).to_radians()))
        .atan()
        .to_degrees();
    println!(
        "   path: level within {PATH_LEVEL:.0} deg of each window, climbing {climb:.1} deg between them"
    );
    // Lying level on its tail: gravity at the centre of mass, buoyancy at its centre, levered
    // about the tail.
    let (tail, _) = b.ends();
    let lift = G * (b.mass * -tail - GLYCEROL_RHO * b.volume * (b.cob.z - tail));
    let mb5 = b.m * 5e-3;
    println!(
        "   standing up: lifting the nose about the tail takes {lift:.2e} N m; at 5 mT a field turning across the swimmer gives at most m B / 2 = {:.2e} averaged over a turn ({:.0}%), a still field m B = {mb5:.2e}",
        0.5 * mb5,
        100.0 * lift / (0.5 * mb5)
    );
    let t_map = Instant::now();
    let map = Map::new(2e-3);
    let coils_f = Field::Coils(&map);
    // The map against cfsem directly, at points through the glycerol.
    {
        let pairs = coils::pairs();
        let mut rng = Rng(7);
        let (mut eb, mut ej) = (0.0_f64, 0.0_f64);
        let (mut bmax, mut jmax) = (0.0_f64, 0.0_f64);
        for _ in 0..200 {
            let p = V3::new(
                (rng.uniform() - 0.5) * 2.0 * LIQ_R * 0.7,
                (rng.uniform() - 0.5) * 2.0 * LIQ_R * 0.7,
                (rng.uniform() - 0.5) * 2.0 * LIQ_H * 0.95,
            );
            let got = map.at(&p);
            let dl = 0.25e-3;
            for k in 0..3 {
                let mut pts = vec![p];
                for ax in 0..3 {
                    let mut e = V3::zeros();
                    e[ax] = dl;
                    pts.push(p + e);
                    pts.push(p - e);
                }
                let bb = pairs[k].field(&pts);
                let mut j = Matrix3::zeros();
                for ax in 0..3 {
                    j.set_column(ax, &((bb[1 + 2 * ax] - bb[2 + 2 * ax]) / (2.0 * dl)));
                }
                eb = eb.max((got[k].0 - bb[0]).norm());
                ej = ej.max((got[k].1 - j).norm());
                bmax = bmax.max(bb[0].norm());
                jmax = jmax.max(j.norm());
            }
        }
        println!(
            "   the coils' field on a 2 mm grid ({:.1} s): against cfsem at 200 points, field within {:.1e} of its size, gradient within {:.1e} of its largest",
            t_map.elapsed().as_secs_f64(),
            eb / bmax,
            ej / jmax
        );
    }
    // The controller: pure pursuit, the aim corrected for sinking, a slow vertical integral; no
    // loop on the seen axis (it made runs worse: below).
    let base = Gains {
        v: 6e-3,
        look: 6e-3,
        kp: 0.0,
        ki: 0.0,
        kz: 0.5,
        comp: true,
        stand: true,
    };
    println!(
        "   controller: {:.0} mm/s along the path, target {:.0} mm ahead, aim corrected for sinking (the model's drag at the seen axis), field rate under {CTRL_KAPPA} of step-out",
        base.v * 1e3,
        base.look * 1e3
    );
    println!(
        "   columns: axis elevation after standing up (deg), when the body left the floor (s), finished (or how far), time, the middle's distance from the path (rms / most), closest to the insert and the glass (mm, negative = pressed in), time touching after lift-off, what it touched, the coils' mean power and the worst pair's heating with no cooling"
    );
    header();
    let mut kept = None;
    for &b0 in &[2e-3, 3e-3, 4e-3, 5e-3] {
        let o = course_run(&b, &coils_f, &c, base, b0, dt, 1, 120.0);
        show(&format!("{:.0} mT, coils", b0 * 1e3), &o);
        if b0 == 5e-3 {
            kept = Some(o);
        }
    }
    let o = course_run(&b, &uni, &c, base, 5e-3, dt, 1, 120.0);
    show("5 mT, uniform field", &o);
    let o = course_run(&b, &coils_f, &c, Gains { kp: 0.5, ki: 0.5, ..base }, 5e-3, dt, 1, 120.0);
    show("5 mT, coils, plus a loop on the axis", &o);
    let o = course_run(&b, &coils_f, &c, Gains { comp: false, ..base }, 5e-3, dt, 1, 120.0);
    show("5 mT, coils, point at the target", &o);
    let o = course_run(&b, &coils_f, &c, Gains { stand: false, ..base }, 5e-3, dt, 1, 120.0);
    show("5 mT, coils, no stand-up", &o);
    let o = course_run(&b, &coils_f, &c, base, 5e-3, 0.25e-3, 1, 120.0);
    show("5 mT, coils, 0.25 ms steps", &o);
    for seed in 2..=5 {
        let o = course_run(&b, &coils_f, &c, base, 5e-3, dt, seed, 120.0);
        show(&format!("5 mT, coils, camera noise seed {seed}"), &o);
    }
    for seed in 2..=3 {
        let o = course_run(&b, &coils_f, &c, base, 4e-3, dt, seed, 120.0);
        show(&format!("4 mT, coils, camera noise seed {seed}"), &o);
    }
    if let Some(o) = kept {
        println!("   the 5 mT run, every 2 s: t, the swimmer's middle (mm), axis elevation (deg), field rate (Hz)");
        let mut next = 0.0;
        for (t, x, a, f) in &o.trace {
            if *t >= next {
                println!(
                    "   {:>6.1}  ({:>6.2}, {:>6.2}, {:>6.2})  {:>7.1}  {:>6.2}",
                    t,
                    x.x * 1e3,
                    x.y * 1e3,
                    x.z * 1e3,
                    a.z.asin().to_degrees(),
                    f
                );
                next += 2.0;
            }
        }
    }
    println!("   ({:.0} s in all)", t_all.elapsed().as_secs_f64());
}

/// `steer wall`: the weightless swimmer swims level into a box wall (5 mT); the deepest contact
/// over time, for the semi-implicit drag at two steps and the plain drag at 20 us.
pub fn run_wall() {
    let b = Body::new(shape::chosen());
    let uni = Field::Uniform;
    println!("steer wall: {} contact spheres", b.spheres.len());
    for (label, dt, explicit, solref) in [
        ("semi-implicit 1 ms, solref 2 steps", 1e-3, false, 2e-3),
        ("semi-implicit 1 ms, solref 20 ms", 1e-3, false, 2e-2),
        ("semi-implicit 0.25 ms, solref 2 steps", 0.25e-3, false, 0.5e-3),
        ("semi-implicit 0.25 ms, solref 20 ms", 0.25e-3, false, 2e-2),
        ("plain 20 us, solref 20 ms", 2e-5, true, 2e-2),
        ("plain 20 us, solref 1 ms", 2e-5, true, 1e-3),
    ] {
        let env = format!(
            "    <geom type=\"box\" pos=\"0.03 0 0\" size=\"0.005 0.02 0.02\" contype=\"0\" conaffinity=\"1\" friction=\"{FRICTION} 0.005 0.0001\" solref=\"{solref} 1\"/>\n"
        );
        let mut s = Sim::new(&b, &uni, dt, &env, V3::zeros(), V3::x(), V3::z(), FRICTION, false, solref);
        s.explicit = explicit;
        s.b0 = 5e-3;
        s.f = KAPPA * b.fso_t_free * 5e-3;
        let mut line = String::new();
        let mut deepest = 0.0_f64;
        let mut vmax = 0.0_f64;
        let mut next = 0.5;
        while s.t < 6.0 {
            s.step();
            let d = s.contacts().iter().map(|k| k.depth).fold(0.0_f64, f64::max);
            deepest = deepest.max(d);
            vmax = vmax.max(s.vel().0.norm());
            if s.t >= next {
                line += &format!(" {:.0}:{:.3}/{:.1}", s.t * 10.0, deepest * 1e3, vmax * 1e3);
                deepest = 0.0;
                vmax = 0.0;
                next += 0.5;
            }
        }
        println!("{label:<40} (t*10: deepest mm / fastest mm/s per 0.5 s){line}");
    }
}


/// `steer contact`: sim-core's contact alone, at an everyday scale: a free 1 kg body (capsule
/// or cylinder lying down, or a sphere; 5 cm radius) dropped from 0.3 m onto fixed ground (a
/// box, a plane, a cylinder standing as a disc, or a cylinder lying on its side, its top at 0),
/// the ground listed before or after the body. 2 ms steps, default contact settings. Resting,
/// the body's centre is at +0.050 m.
pub fn run_contact() {
    let ground = |kind: &str| match kind {
        "box" => "<geom type=\"box\" pos=\"0 0 -0.1\" size=\"1 1 0.1\"/>",
        "plane" => "<geom type=\"plane\" size=\"1 1 0.01\"/>",
        "disc" => "<geom type=\"cylinder\" pos=\"0 0 -0.1\" size=\"1 0.1\"/>",
        _ => "<geom type=\"cylinder\" pos=\"0 0 -0.1\" euler=\"90 0 0\" size=\"0.1 1\"/>",
    };
    println!("steer contact: a free 1 kg body dropped onto fixed ground (resting centre +0.050 m)");
    println!(
        "   {:<9} {:<14} {:>12} {:>12} {:>10} {:>10}",
        "body", "ground", "listed", "centre m", "contacts", "deepest m"
    );
    for shape in ["sphere", "capsule", "cylinder"] {
        for kind in ["plane", "box", "disc", "cylinder side"] {
            for ground_first in [true, false] {
                let body = format!(
                    "<body pos=\"0 0 0.3\" euler=\"0 90 0\"><freejoint/><geom type=\"{shape}\" size=\"0.05 0.1\" mass=\"1\"/></body>"
                );
                let g = ground(kind);
                let world = if ground_first {
                    format!("{g}\n{body}")
                } else {
                    format!("{body}\n<body name=\"gb\">{g}</body>")
                };
                let xml = format!(
                    "<mujoco><option timestep=\"0.002\"/><worldbody>\n{world}\n</worldbody></mujoco>"
                );
                let model = load_model(&xml).expect("contact model");
                let mut d = model.make_data();
                for _ in 0..1500 {
                    d.step(&model).expect("step");
                }
                let depth = d.contacts[..d.ncon].iter().map(|k| k.depth).fold(0.0_f64, f64::max);
                let ok = (d.qpos[2] - 0.05).abs() < 2e-3;
                println!(
                    "   {shape:<9} {kind:<14} {:>12} {:>+12.4} {:>10} {:>10.4}  {}",
                    if ground_first { "ground first" } else { "body first" },
                    d.qpos[2],
                    d.ncon,
                    depth,
                    if ok { "rests" } else { "FAILS" }
                );
            }
        }
    }
}
