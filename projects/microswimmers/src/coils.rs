//! `coils`: three nested Helmholtz pairs round the jar (coil-designer's turn-by-turn model): do
//! they fit, what driving them takes, how round and even the rotating field is over the
//! glycerol, and how hard the field's unevenness pulls the swimmer's magnet against its swimming.
//!
//! A field turning in the plane normal to n (the way the swimmer should go) is
//! B(t) = B0 (e1 cos wt + e2 sin wt) at the centre, from pair currents I_k = B0 (e1_k cos wt +
//! e2_k sin wt) / k_k (k_k the pair's centre field per amp). Elsewhere it is P cos wt + Q sin wt:
//! an ellipse, tilted. The magnet trails the field by phi, sin phi = (f / f_so) B0 / |B|, and
//! feels F = (m . grad) B, averaged over a turn.

use std::f64::consts::PI;

use coil_designer::coil::{Coil, Pack, Wire};
use nalgebra::{Matrix3, SymmetricEigen};
use rayon::prelude::*;

use crate::check::GLYCEROL_MU;
use crate::helix::V3;
use crate::rft::swim_locked;
use crate::swimmer::Swimmer;

/// (J/(kg K)) copper.
const C_CU: f64 = 385.0;
/// The jar, upright on z: glass outer radius and half height, and the glycerol inside.
const JAR_R: f64 = 28e-3;
const JAR_H: f64 = 30e-3;
const LIQ_R: f64 = 25e-3;
const LIQ_H: f64 = 25e-3;
/// Radius of the central region reported on its own.
const CORE: f64 = 10e-3;
/// Operate at this fraction of step-out when comparing the pull with swimming.
const KAPPA: f64 = 0.75;

struct Pair {
    name: &'static str,
    axis: V3,
    u: V3,
    v: V3,
    coil: Coil,
    r_mean: f64,
    spacing: f64,
    turns: usize,
}

impl Pair {
    fn new(name: &'static str, axis: V3, u: V3, r_mean: f64, wire: Wire, n: usize) -> Pair {
        let depth = n as f64 * wire.outer;
        let pack = |z: f64| Pack {
            r_in: r_mean - 0.5 * depth,
            z_mid: z,
            layers: n,
            per_layer: n,
            wire,
            sign: 1.0,
            partial: 0,
        };
        Pair {
            name,
            axis,
            u,
            v: axis.cross(&u),
            coil: Coil {
                packs: vec![pack(-0.5 * r_mean), pack(0.5 * r_mean)],
            },
            r_mean,
            spacing: r_mean,
            turns: n * n,
        }
    }

    /// Field per amp at world points.
    fn field(&self, pts: &[V3]) -> Vec<V3> {
        let mut r = Vec::with_capacity(pts.len());
        let mut z = Vec::with_capacity(pts.len());
        let mut er = Vec::with_capacity(pts.len());
        for p in pts {
            let zl = p.dot(&self.axis);
            let rv = p - self.axis * zl;
            let rl = rv.norm();
            r.push(rl);
            z.push(zl);
            er.push(if rl > 0.0 { rv / rl } else { V3::zeros() });
        }
        let (br, bz) = self.coil.field(1.0, &r, &z);
        (0..pts.len())
            .map(|i| er[i] * br[i] + self.axis * bz[i])
            .collect()
    }

    /// Points on the winding packs' outer surfaces (each pack's section perimeter, swept round).
    fn surface(&self) -> Vec<V3> {
        let mut out = Vec::new();
        for p in &self.coil.packs {
            let (r0, r1) = (p.r_in, p.r_in + p.depth());
            let (z0, z1) = (p.z_mid - 0.5 * p.width(), p.z_mid + 0.5 * p.width());
            let k = 8;
            let mut rect = Vec::new();
            for i in 0..k {
                let s = i as f64 / k as f64;
                rect.push((r0 + s * (r1 - r0), z0));
                rect.push((r1, z0 + s * (z1 - z0)));
                rect.push((r1 - s * (r1 - r0), z1));
                rect.push((r0, z1 - s * (z1 - z0)));
            }
            for j in 0..360 {
                let a = 2.0 * PI * j as f64 / 360.0;
                let dir = self.u * a.cos() + self.v * a.sin();
                for &(r, z) in &rect {
                    out.push(dir * r + self.axis * z);
                }
            }
        }
        out
    }
}

fn min_distance(a: &[V3], b: &[V3]) -> f64 {
    a.par_iter()
        .map(|p| b.iter().map(|q| (p - q).norm()).fold(f64::INFINITY, f64::min))
        .reduce(|| f64::INFINITY, f64::min)
}

/// Distance from a point to the solid jar (negative inside).
fn jar_distance(p: &V3) -> f64 {
    let rho = (p.x * p.x + p.y * p.y).sqrt();
    let (dr, dz) = (rho - JAR_R, p.z.abs() - JAR_H);
    if dr <= 0.0 && dz <= 0.0 {
        dr.max(dz)
    } else {
        (dr.max(0.0).powi(2) + dz.max(0.0).powi(2)).sqrt()
    }
}

/// Room left for the jar lifted straight up and out (or lowered in): the winding points above
/// the jar's bottom, measured sideways from the jar's axis, less the jar's radius.
fn lift_clearance(p: &V3) -> f64 {
    if p.z < -JAR_H {
        f64::INFINITY
    } else {
        (p.x * p.x + p.y * p.y).sqrt() - JAR_R
    }
}

/// e1, e2 spanning the plane normal to n (right-handed: e1 x e2 = n).
fn plane(n: &V3) -> (V3, V3) {
    let t = if n.x.abs() < 0.9 { V3::x() } else { V3::y() };
    let e1 = (t - n * n.dot(&t)).normalize();
    (e1, n.cross(&e1))
}

pub fn run() {
    let wire = Wire::awg(20);
    let n = 13;
    let pairs = [
        Pair::new("x", V3::x(), V3::y(), 75e-3, wire, n),
        Pair::new("y", V3::y(), V3::z(), 92e-3, wire, n),
        Pair::new("z", V3::z(), V3::x(), 112e-3, wire, n),
    ];
    let k: Vec<f64> = pairs
        .iter()
        .map(|p| p.field(&[V3::zeros()])[0].dot(&p.axis))
        .collect();

    println!(
        "coils: three nested Helmholtz pairs round a {:.0} mm jar, {:.0} mm tall (glycerol {:.0} x {:.0} mm)",
        2.0 * JAR_R * 1e3,
        2.0 * JAR_H * 1e3,
        2.0 * LIQ_R * 1e3,
        2.0 * LIQ_H * 1e3
    );
    println!(
        "1. The coils: each coil AWG {} ({:.3} mm bare), {n} x {n} = {} turns, pack {:.1} x {:.1} mm; spacing = mean radius",
        wire.awg,
        wire.bare * 1e3,
        n * n,
        n as f64 * wire.outer * 1e3,
        n as f64 * wire.outer * 1e3
    );
    println!(
        "   {:>4} {:>8} {:>9} {:>10} {:>8} {:>8} {:>8} {:>8} {:>9}",
        "pair", "R mm", "mT/A", "ohm", "mH", "L/R ms", "wire m", "Cu kg", "turns"
    );
    let mut res = Vec::new();
    let mut ind = Vec::new();
    let mut cu = Vec::new();
    for (p, &kk) in pairs.iter().zip(&k) {
        let (r, l, m) = (
            p.coil.resistance(),
            p.coil.inductance(),
            p.coil.copper_mass(),
        );
        println!(
            "   {:>4} {:>8.1} {:>9.4} {:>10.3} {:>8.3} {:>8.3} {:>8.1} {:>8.3} {:>9}",
            p.name,
            p.r_mean * 1e3,
            kk * 1e3,
            r,
            l * 1e3,
            l / r * 1e3,
            p.coil.wire_length(),
            m,
            2 * p.turns
        );
        res.push(r);
        ind.push(l);
        cu.push(m);
    }
    let surf: Vec<Vec<V3>> = pairs.iter().map(|p| p.surface()).collect();
    for (i, j) in [(0, 1), (1, 2), (0, 2)] {
        println!(
            "   clearance between the {} and {} windings: {:.1} mm",
            pairs[i].name,
            pairs[j].name,
            min_distance(&surf[i], &surf[j]) * 1e3
        );
    }
    for (p, s) in pairs.iter().zip(&surf) {
        let d = s
            .iter()
            .map(jar_distance)
            .fold(f64::INFINITY, f64::min);
        let lift = s
            .iter()
            .map(lift_clearance)
            .fold(f64::INFINITY, f64::min);
        println!(
            "   the {} windings: {:.1} mm from the jar where it stands, {:.1} mm to spare lifting it out (spacing {:.1} mm)",
            p.name,
            d * 1e3,
            lift * 1e3,
            p.spacing * 1e3
        );
    }
    println!(
        "   (orthogonal pairs on one centre have no mutual inductance: each pair's field through the other's loops cancels by symmetry)"
    );

    // The swimmer: step-out per tesla, speed per Hz, translation block.
    let sw = Swimmer::desk();
    let model = sw.model(GLYCEROL_MU);
    let held = swim_locked(&model.r);
    let m = sw.moment();
    let fso_t = m / (2.0 * PI * held.torque); // Hz per tesla
    let v_hz = held.u.z * 2.0 * PI;
    let a_body: Matrix3<f64> = model.r.fixed_view::<3, 3>(0, 0).into_owned();
    println!(
        "   swimmer (from swim): step-out {:.3} Hz per mT, {:.3} mm/s per Hz",
        fso_t * 1e-3,
        v_hz * 1e3
    );

    // -----------------------------------------------------------------------------------------
    let planes = [
        ("climb (n = z)", V3::z()),
        ("level along x", V3::x()),
        ("level along y", V3::y()),
    ];
    println!("2. Driving it: peak current and voltage (at step-out) per pair, mean power, heating with no cooling");
    println!(
        "   {:>4} {:>14} {:>7} {:>22} {:>22} {:>8} {:>13} {:>10}",
        "B mT", "swim", "f_so", "peak A (x, y, z)", "peak V (x, y, z)", "mean W", "worst C/min", "A/mm^2"
    );
    for b0 in [2e-3, 5e-3] {
        let f = fso_t * b0;
        for (label, nv) in &planes {
            let (e1, e2) = plane(nv);
            let mut ia = [0.0; 3];
            let mut va = [0.0; 3];
            let mut watts = 0.0;
            let mut heat = 0.0_f64;
            let mut jmax = 0.0_f64;
            for kk in 0..3 {
                let amp = (e1[kk].powi(2) + e2[kk].powi(2)).sqrt();
                let i_pk = b0 * amp / k[kk];
                let w = 0.5 * res[kk] * i_pk * i_pk;
                ia[kk] = i_pk;
                va[kk] = i_pk * (res[kk].powi(2) + (2.0 * PI * f * ind[kk]).powi(2)).sqrt();
                watts += w;
                heat = heat.max(w / (cu[kk] * C_CU) * 60.0);
                jmax = jmax.max(i_pk / (PI * wire.bare.powi(2) / 4.0) * 1e-6);
            }
            println!(
                "   {:>4.0} {label:>14} {f:>7.2} {:>22} {:>22} {watts:>8.1} {heat:>13.2} {jmax:>10.2}",
                b0 * 1e3,
                format!("{:.2}, {:.2}, {:.2}", ia[0], ia[1], ia[2]),
                format!("{:.1}, {:.1}, {:.1}", va[0], va[1], va[2]),
            );
        }
    }

    // -----------------------------------------------------------------------------------------
    // Sample points in the glycerol, with neighbours for finite differences.
    let g = 13;
    let mut pts = Vec::new();
    for i in 0..g {
        for j in 0..g {
            for l in 0..g {
                let x = -LIQ_R + 2.0 * LIQ_R * i as f64 / (g - 1) as f64;
                let y = -LIQ_R + 2.0 * LIQ_R * j as f64 / (g - 1) as f64;
                let z = -LIQ_H + 2.0 * LIQ_H * l as f64 / (g - 1) as f64;
                if x * x + y * y <= LIQ_R * LIQ_R * (1.0 + 1e-9) {
                    pts.push(V3::new(x, y, z));
                }
            }
        }
    }
    let dl = 0.25e-3;
    let mut all = Vec::with_capacity(pts.len() * 7);
    for p in &pts {
        all.push(*p);
        for ax in 0..3 {
            let mut e = V3::zeros();
            e[ax] = dl;
            all.push(p + e);
            all.push(p - e);
        }
    }
    let b: Vec<Vec<V3>> = pairs.iter().map(|p| p.field(&all)).collect();

    println!(
        "3. The rotating field over the glycerol ({} points; core = within {:.0} mm of the centre), per unit B0",
        pts.len(),
        CORE * 1e3
    );
    println!(
        "   {:>14} {:>6} {:>17} {:>13} {:>13} {:>10}",
        "swim", "where", "mean |B| min-max", "weakest |B|", "ellipticity", "tilt deg"
    );
    let mut pull_rows = Vec::new();
    for (label, nv) in &planes {
        let (e1, e2) = plane(nv);
        // P and Q at sample index s.
        let pq = |s: usize| -> (V3, V3) {
            let mut p = V3::zeros();
            let mut q = V3::zeros();
            for kk in 0..3 {
                p += b[kk][s] * (e1[kk] / k[kk]);
                q += b[kk][s] * (e2[kk] / k[kk]);
            }
            (p, q)
        };
        let rot = Matrix3::from_columns(&[e1, e2, *nv]);
        let a_world = rot * a_body * rot.transpose();
        let a_inv = a_world.try_inverse().expect("drag invertible");
        let swim_speed = v_hz * KAPPA * fso_t; // per tesla of B0
        let nt = 64;
        for (where_, core_only) in [("all", false), ("core", true)] {
            let (mut mean_lo, mut mean_hi) = (f64::INFINITY, 0.0_f64);
            let (mut weakest, mut ellip, mut tilt) = (f64::INFINITY, 0.0_f64, 0.0_f64);
            let mut pull = 0.0_f64;
            for (ip, p) in pts.iter().enumerate() {
                if core_only && p.norm() > CORE {
                    continue;
                }
                let s0 = 7 * ip;
                let (p0, q0) = pq(s0);
                let gram = Matrix3::new(
                    p0.dot(&p0), p0.dot(&q0), 0.0,
                    p0.dot(&q0), q0.dot(&q0), 0.0,
                    0.0, 0.0, 0.0,
                );
                let ev = SymmetricEigen::new(gram).eigenvalues;
                let mut evs: Vec<f64> = ev.iter().copied().filter(|v| *v > 1e-30).collect();
                evs.sort_by(f64::total_cmp);
                let (lo, hi) = (evs[0].sqrt(), evs[evs.len() - 1].sqrt());
                weakest = weakest.min(lo);
                ellip = ellip.max(hi / lo - 1.0);
                let axis = p0.cross(&q0).normalize();
                tilt = tilt.max(axis.dot(nv).clamp(-1.0, 1.0).acos().to_degrees());
                // Jacobians of P and Q by central differences: column a = d/dx_a.
                let mut jp = Matrix3::zeros();
                let mut jq = Matrix3::zeros();
                for ax in 0..3 {
                    let (pp, qp) = pq(s0 + 1 + 2 * ax);
                    let (pm, qm) = pq(s0 + 2 + 2 * ax);
                    jp.set_column(ax, &((pp - pm) / (2.0 * dl)));
                    jq.set_column(ax, &((qp - qm) / (2.0 * dl)));
                }
                let mut mean = 0.0;
                let mut force = V3::zeros();
                for it in 0..nt {
                    let wt = 2.0 * PI * it as f64 / nt as f64;
                    let bt = p0 * wt.cos() + q0 * wt.sin();
                    let bn = bt.norm();
                    mean += bn / nt as f64;
                    // The magnet trails the field by phi about the field's turning axis.
                    let sphi = (KAPPA / bn).min(1.0);
                    let bh = bt / bn;
                    let mh = bh * (1.0 - sphi * sphi).sqrt() - axis.cross(&bh) * sphi;
                    let jt = jp * wt.cos() + jq * wt.sin();
                    // (m . grad) B, with grad B = J (J_ij = dB_i/dx_j): F_i = sum_j m_j dB_i/dx_j
                    // = sum_j m_j dB_j/dx_i for a curl-free field; use J^T m.
                    force += jt.transpose() * mh / nt as f64;
                }
                mean_lo = mean_lo.min(mean);
                mean_hi = mean_hi.max(mean);
                // Per tesla of B0: force m B0 * force_hat; drift = A^-1 F.
                let drift = (a_inv * force * m).norm();
                pull = pull.max(drift / swim_speed);
            }
            println!(
                "   {label:>14} {where_:>6} {:>17} {weakest:>13.3} {ellip:>13.3} {tilt:>10.2}",
                format!("{mean_lo:.3}-{mean_hi:.3}"),
            );
            pull_rows.push((label.to_string(), where_, pull));
        }
    }
    println!("4. The magnet's pull from the unevenness, against swimming, at {KAPPA} of step-out");
    for (label, where_, pull) in &pull_rows {
        println!("   {label:>14} {where_:>6}: worst drift / swimming speed {pull:.3}");
    }
}
