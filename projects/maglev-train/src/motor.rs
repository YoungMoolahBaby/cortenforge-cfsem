//! `motor`: a long-stator linear synchronous motor that pushes on the cart's own lift magnets.
//!
//! The lift magnets alternate polarity along each row (north at the rear corner, south at the
//! front, 100 mm apart): the Transrapid trick of letting the levitation magnets double as the
//! motor's field. A fine-pitch flat three-phase stator on the rail's underside (a printed-circuit
//! coil board, in the gap), one coil track under each magnet row with phases reversed between
//! rows as the rows' polarities are, pushes on them without contact. A first design with the coil
//! pitch equal to the magnet spacing (100 mm) gave thrust only in narrow spikes, when a 12 mm
//! magnet sat right under a coil side; see TAU.
//!
//! Coils are cfsem linear filaments mirrored in the steel; the magnets are their equivalent
//! current loops. Thrust is checked two ways: the Lorentz force on the magnet loops in the coils'
//! field, and the coenergy route F = I dPsi/dx with Psi the magnets' flux through a coil
//! (`mutual_inductance_circular_to_linear`).

use std::f64::consts::PI;

use cfsem::physics::circular_filament::mutual_inductance_circular_to_linear_scalar;
use cfsem::physics::linear_filament::flux_density_linear_filament;

use crate::em::magnet_loops;

/// (m) stator pole pitch. The two magnets in a row are 100 mm = 7 TAU apart, an odd multiple, so
/// they sit pi apart in electrical angle, which their opposite polarities need; and TAU is close
/// to the magnet size, so some coil side is always under each magnet.
pub const TAU: f64 = 0.1 / 7.0;
pub const LX: f64 = 0.05;
pub const LY: f64 = 0.025;
const CORNERS: [(f64, f64); 4] = [(-LX, -LY), (LX, -LY), (LX, LY), (-LX, LY)];
/// Alternating polarity, as in `stray`: corners 0 and 2 north-up, 1 and 3 south-up.
const POLARITY: [f64; 4] = [1.0, -1.0, 1.0, -1.0];
const COIL_W: f64 = 0.024; // m, across the track
const COIL_Z: f64 = -0.00075; // m, coil plane just under the steel face
pub const TURNS: f64 = 20.0;
/// (ohm) one coil: 20 turns of AWG 24 (84.2 mohm/m) round 2 (TAU + 24 mm).
pub const COIL_R: f64 = 20.0 * 2.0 * (TAU + 0.024) * 0.0842;
const COM_BELOW_FACE: f64 = 0.015;
const LOOPS: usize = 24;

/// One stator coil: a rectangle on the coil plane, its phase (0, 1, 2) and winding sense.
#[derive(Clone, Copy)]
struct StatorCoil {
    x0: f64,
    yc: f64,
    phase: usize,
    sign: f64,
}

/// Coils along the track between x_from and x_to (both rows, all phases).
fn stator(x_from: f64, x_to: f64) -> Vec<StatorCoil> {
    let mut v = Vec::new();
    let mut n = ((x_from - 2.0 * TAU) / (2.0 * TAU)).floor() as i64;
    loop {
        let base = n as f64 * 2.0 * TAU;
        if base > x_to + 2.0 * TAU {
            break;
        }
        for phase in 0..3 {
            let x0 = base + phase as f64 * 2.0 * TAU / 3.0;
            if x0 + TAU < x_from || x0 > x_to {
                continue;
            }
            for (yc, sign) in [(-LY, 1.0), (LY, -1.0)] {
                v.push(StatorCoil {
                    x0,
                    yc,
                    phase,
                    sign,
                });
            }
        }
        n += 1;
    }
    v
}

/// The rectangle's corners, counter-clockwise seen from below (+z up).
fn corners(c: &StatorCoil, z: f64) -> [(f64, f64, f64); 4] {
    let (y0, y1) = (c.yc - 0.5 * COIL_W, c.yc + 0.5 * COIL_W);
    [
        (c.x0, y0, z),
        (c.x0 + TAU, y0, z),
        (c.x0 + TAU, y1, z),
        (c.x0, y1, z),
    ]
}

/// Straight segments (start, delta, ampere-turns) of the coils and their images (z -> -z, same
/// current: iron), for phase currents `i`.
/// A straight current segment: start, delta (m), ampere-turns.
type Segment = ((f64, f64, f64), (f64, f64, f64), f64);

fn segments(coils: &[StatorCoil], i: [f64; 3]) -> Vec<Segment> {
    let mut v = Vec::new();
    for c in coils {
        let cur = c.sign * TURNS * i[c.phase];
        if cur == 0.0 {
            continue;
        }
        for z in [COIL_Z, -COIL_Z] {
            let p = corners(c, z);
            for k in 0..4 {
                let (a, b) = (p[k], p[(k + 1) % 4]);
                v.push((a, (b.0 - a.0, b.1 - a.1, b.2 - a.2), cur));
            }
        }
    }
    v
}

/// Force (x, z) and moments (about x, about y) on the cart, its centre of mass COM_BELOW_FACE
/// below the magnet faces, from the coils carrying phase currents `i`, by the Lorentz force on the
/// magnets' loops (each split into `seg` points).
pub fn wrench(cart_x: f64, gap: f64, i: [f64; 3]) -> [f64; 4] {
    let coils = stator(cart_x - 0.35, cart_x + 0.35);
    let segs = segments(&coils, i);
    let (sx, sy, sz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
        segs.iter().map(|s| s.0.0).collect(),
        segs.iter().map(|s| s.0.1).collect(),
        segs.iter().map(|s| s.0.2).collect(),
    );
    let (dx, dy, dz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
        segs.iter().map(|s| s.1.0).collect(),
        segs.iter().map(|s| s.1.1).collect(),
        segs.iter().map(|s| s.1.2).collect(),
    );
    let cur: Vec<f64> = segs.iter().map(|s| s.2).collect();
    let zero = vec![0.0; segs.len()];
    let seg = 64;
    let com = (cart_x, 0.0, -gap - COM_BELOW_FACE);
    let mut w = [0.0; 4];
    for (k, &(cx, cy)) in CORNERS.iter().enumerate() {
        let (loops, il) = magnet_loops(0.006, 0.004, 1.3, LOOPS, gap);
        let il = il * POLARITY[k];
        for &(r, z) in &loops {
            let pts: Vec<(f64, f64, f64)> = (0..seg)
                .map(|j| {
                    let t = 2.0 * PI * (j as f64 + 0.5) / seg as f64;
                    (cart_x + cx + r * t.cos(), cy + r * t.sin(), z)
                })
                .collect();
            let (px, py, pz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
                pts.iter().map(|p| p.0).collect(),
                pts.iter().map(|p| p.1).collect(),
                pts.iter().map(|p| p.2).collect(),
            );
            let (mut bx, mut by, mut bz) = (vec![0.0; seg], vec![0.0; seg], vec![0.0; seg]);
            flux_density_linear_filament(
                (&px, &py, &pz),
                (&sx, &sy, &sz),
                (&dx, &dy, &dz),
                &cur,
                &zero,
                (&mut bx, &mut by, &mut bz),
            )
            .expect("lengths");
            for j in 0..seg {
                let t = 2.0 * PI * (j as f64 + 0.5) / seg as f64;
                let dl = (
                    -t.sin() * 2.0 * PI * r / seg as f64,
                    t.cos() * 2.0 * PI * r / seg as f64,
                    0.0,
                );
                // dF = I dl x B
                let f = (
                    il * (dl.1 * bz[j] - dl.2 * by[j]),
                    il * (dl.2 * bx[j] - dl.0 * bz[j]),
                    il * (dl.0 * by[j] - dl.1 * bx[j]),
                );
                let rr = (pts[j].0 - com.0, pts[j].1 - com.1, pts[j].2 - com.2);
                w[0] += f.0;
                w[1] += f.2;
                w[2] += rr.1 * f.2 - rr.2 * f.1; // moment about x
                w[3] += rr.2 * f.0 - rr.0 * f.2; // moment about y
            }
        }
    }
    w
}

/// Flux linkage (Wb) of one phase's coils (unit current sense, with turns) from the magnets and
/// their images, by `mutual_inductance_circular_to_linear` on coil sides split into short pieces
/// (its docs: accurate only for segments much shorter than the loop radius).
fn linkage(cart_x: f64, gap: f64, phase: usize, piece: f64, true_mid: bool) -> f64 {
    let coils: Vec<StatorCoil> = stator(cart_x - 0.35, cart_x + 0.35)
        .into_iter()
        .filter(|c| c.phase == phase)
        .collect();
    let mut psi = 0.0;
    for (k, &(cx, cy)) in CORNERS.iter().enumerate() {
        let (loops, il) = magnet_loops(0.006, 0.004, 1.3, LOOPS, gap);
        let il = il * POLARITY[k];
        let (mx, my) = (cart_x + cx, cy);
        for c in &coils {
            let p = corners(c, COIL_Z);
            for e in 0..4 {
                let (a, b) = (p[e], p[(e + 1) % 4]);
                let len = (b.0 - a.0).hypot(b.1 - a.1);
                let n = (len / piece).ceil() as usize;
                for s in 0..n {
                    let (t0, t1) = (s as f64 / n as f64, (s + 1) as f64 / n as f64);
                    let p0 = (
                        a.0 + t0 * (b.0 - a.0) - mx,
                        a.1 + t0 * (b.1 - a.1) - my,
                        a.2,
                    );
                    let p1 = (
                        a.0 + t1 * (b.0 - a.0) - mx,
                        a.1 + t1 * (b.1 - a.1) - my,
                        a.2,
                    );
                    for &(r, z) in &loops {
                        // The magnet loop and its image (z -> -z, same current).
                        for zz in [z, -z] {
                            let m = if true_mid {
                                mutual_midpoint((r, zz, 1.0), p0, p1)
                            } else {
                                mutual_inductance_circular_to_linear_scalar((r, zz, 1.0), p0, p1)
                            };
                            psi += c.sign * TURNS * il * m;
                        }
                    }
                }
            }
        }
    }
    psi
}

/// Commutation: phase currents for thrust current `iq` at cart position x, sinusoidal in the
/// electrical angle pi x / tau with phases 2 pi / 3 apart and the offset `delta` that puts the
/// current in step with the magnets.
pub fn phase_currents(iq: f64, x: f64, delta: f64) -> [f64; 3] {
    let th = PI * x / TAU + delta;
    [0, 1, 2].map(|k| iq * (th - 2.0 * PI * k as f64 / 3.0).cos())
}

/// Thrust per amp of thrust current against cart position over one electrical period, and the
/// other wrench components, at gap g: a table for the cart simulation.
pub struct MotorTable {
    pub xs: Vec<f64>,
    pub w: Vec<[f64; 4]>,
    pub delta: f64,
    pub kf: f64,
    /// (W per A^2 of thrust current) mean stator copper loss, with block switching: the coils
    /// under the cart (its 140 mm length plus 10 mm each end) carry the commutated currents.
    pub copper: f64,
}

impl MotorTable {
    pub fn build(gap: f64, n: usize) -> Self {
        // Phase A alone: thrust F_A(x) ~ K sin(pi x / tau + phi). Fit phi to commutate.
        let samples: Vec<(f64, f64)> = (0..24)
            .map(|k| {
                let x = 2.0 * TAU * k as f64 / 24.0;
                (x, wrench(x, gap, [1.0, 0.0, 0.0])[0])
            })
            .collect();
        let (mut s, mut c) = (0.0, 0.0);
        for &(x, f) in &samples {
            s += f * (PI * x / TAU).sin();
            c += f * (PI * x / TAU).cos();
        }
        // F_A = A sin(th) + B cos(th) = R cos(th - atan2(A, B)); current i_A = cos(th + delta)
        // in phase with F_A means delta = -atan2(A, B).
        let delta = -s.atan2(c);
        let xs: Vec<f64> = (0..n).map(|k| 2.0 * TAU * k as f64 / n as f64).collect();
        let w: Vec<[f64; 4]> = xs
            .iter()
            .map(|&x| wrench(x, gap, phase_currents(1.0, x, delta)))
            .collect();
        let kf = w.iter().map(|v| v[0]).sum::<f64>() / n as f64;
        let copper = xs
            .iter()
            .map(|&x| {
                let i = phase_currents(1.0, x, delta);
                stator(x - 0.08, x + 0.08)
                    .iter()
                    .map(|c| i[c.phase] * i[c.phase] * COIL_R)
                    .sum::<f64>()
            })
            .sum::<f64>()
            / xs.len() as f64;
        MotorTable {
            xs,
            w,
            delta,
            kf,
            copper,
        }
    }
    /// Wrench per amp at cart position x (periodic, linear interpolation).
    pub fn at(&self, x: f64) -> [f64; 4] {
        let p = 2.0 * TAU;
        let u = (x.rem_euclid(p)) / p * self.xs.len() as f64;
        let i = u.floor() as usize % self.xs.len();
        let j = (i + 1) % self.xs.len();
        let t = u - u.floor();
        [0, 1, 2, 3].map(|c| self.w[i][c] * (1.0 - t) + self.w[j][c] * t)
    }
}

/// The same mutual inductance with the vector potential taken at the segment's true (Cartesian)
/// midpoint, for comparison with cfsem's `mutual_inductance_circular_to_linear_scalar`.
fn mutual_midpoint(rzn: (f64, f64, f64), p0: (f64, f64, f64), p1: (f64, f64, f64)) -> f64 {
    let m = (
        0.5 * (p0.0 + p1.0),
        0.5 * (p0.1 + p1.1),
        0.5 * (p0.2 + p1.2),
    );
    let (r, phi) = (m.0.hypot(m.1), m.1.atan2(m.0));
    let a =
        cfsem::physics::circular_filament::vector_potential_circular_filament_scalar(rzn, (r, m.2));
    a * (-phi.sin() * (p1.0 - p0.0) + phi.cos() * (p1.1 - p0.1))
}

/// Debug probe: per-phase thrust and lift at 1 A against cart position.
fn phase_curves(gap: f64) {
    println!("   x mm   F_A      F_B      F_C   (thrust, mN/A)    lift A / B / C (mN/A)");
    for k in 0..12 {
        let x = 2.0 * TAU * k as f64 / 12.0;
        let w: Vec<[f64; 4]> = (0..3)
            .map(|p| {
                let mut i = [0.0; 3];
                i[p] = 1.0;
                wrench(x, gap, i)
            })
            .collect();
        println!(
            "   {:>4.0}  {:+7.2}  {:+7.2}  {:+7.2}        {:+7.1} / {:+7.1} / {:+7.1}",
            x * 1e3,
            w[0][0] * 1e3,
            w[1][0] * 1e3,
            w[2][0] * 1e3,
            w[0][1] * 1e3,
            w[1][1] * 1e3,
            w[2][1] * 1e3
        );
    }
}

pub fn run() -> MotorTable {
    let gap = 0.005964;
    println!(
        "0. cfsem's mutual_inductance_circular_to_linear_scalar vs the same formula at the true midpoint"
    );
    let loop_ = (0.01, 0.0, 1.0);
    for (what, p0, p1) in [
        (
            "tangential, r 20 mm, 10 mm long",
            (0.02, -0.005, -0.003),
            (0.02, 0.005, -0.003),
        ),
        (
            "diagonal outward, 10 mm long",
            (0.02, -0.005, -0.003),
            (0.027, 0.002, -0.003),
        ),
        (
            "diagonal inward, 10 mm long",
            (0.027, 0.002, -0.003),
            (0.02, -0.005, -0.003),
        ),
        (
            "crossing the -x axis, 0.2 mm long",
            (-0.02, 0.0001, -0.003),
            (-0.02, -0.0001, -0.003),
        ),
    ] {
        let c = mutual_inductance_circular_to_linear_scalar(loop_, p0, p1);
        let m = mutual_midpoint(loop_, p0, p1);
        println!(
            "   {what:<36} cfsem {c:+.6e} H, true midpoint {m:+.6e} H ({:+.2e})",
            c / m - 1.0
        );
    }
    phase_curves(gap);
    println!(
        "Stator: full-pitch flat coils {:.0} x {:.0} mm, {TURNS:.0} turns ({:.2} ohm), 3 phases 2 tau / 3 apart, \
         one track under each magnet row, {:.2} mm under the steel; cart at gap {:.3} mm",
        TAU * 1e3,
        COIL_W * 1e3,
        COIL_R,
        -COIL_Z * 1e3,
        gap * 1e3
    );
    println!(
        "1. Thrust from phase A at 1 A: Lorentz on the magnet loops vs coenergy I dPsi/dx (circular-to-linear mutual)"
    );
    println!("   (coenergy with cfsem's segment midpoint, and with the true midpoint)");
    for x in [0.002, 0.005, 0.011, 0.019] {
        let f = wrench(x, gap, [1.0, 0.0, 0.0])[0];
        let h = 1e-5;
        for piece in [1e-3, 2.5e-4] {
            let d = |mid: bool| {
                (linkage(x + h, gap, 0, piece, mid) - linkage(x - h, gap, 0, piece, mid))
                    / (2.0 * h)
            };
            let (ec, et) = (d(false), d(true));
            println!(
                "   x {:>4.1} mm, pieces {:.2} mm: Lorentz {f:+.5e} N; coenergy cfsem {ec:+.5e} ({:+.1e}), true midpoint {et:+.5e} ({:+.1e})",
                x * 1e3,
                piece * 1e3,
                ec / f - 1.0,
                et / f - 1.0
            );
        }
    }
    let table = MotorTable::build(gap, 48);
    let fx: Vec<f64> = table.w.iter().map(|v| v[0]).collect();
    let (lo, hi) = fx
        .iter()
        .fold((f64::MAX, f64::MIN), |a, &v| (a.0.min(v), a.1.max(v)));
    let fz: Vec<f64> = table.w.iter().map(|v| v[1]).collect();
    let my: Vec<f64> = table.w.iter().map(|v| v[3]).collect();
    let span = |v: &[f64]| {
        v.iter()
            .fold((f64::MAX, f64::MIN), |a, &x| (a.0.min(x), a.1.max(x)))
    };
    println!(
        "2. Commutated three-phase thrust (sinusoidal in pi x / tau, offset {:.3} rad), per amp of thrust current",
        table.delta
    );
    println!(
        "   thrust constant {:.4} N/A, ripple {:.2}% peak to peak; lift {:+.2e} .. {:+.2e} N/A; pitch moment {:+.2e} .. {:+.2e} N m/A",
        table.kf,
        (hi - lo) / table.kf * 100.0,
        span(&fz).0,
        span(&fz).1,
        span(&my).0,
        span(&my).1
    );
    // Copper: two coils per phase per period energized round the cart (block switching).
    let p = table.copper;
    println!(
        "   copper loss at 1 A, coils under the cart energized: {p:.2} W; \
         thrust per watt {:.3} N/W",
        table.kf / p
    );
    table
}
