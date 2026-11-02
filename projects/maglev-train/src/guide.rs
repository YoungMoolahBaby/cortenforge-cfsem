//! `guide`: contactless sideways guidance, a second coil layer on the stator board.
//!
//! The steel plane pulls straight up only, so sideways force has to come from track currents
//! pushing on the lift magnets. Guidance coils span both magnet rows: 50 mm across, so their long
//! sides run along the track right over the rows. The rows have opposite polarity and those sides
//! carry opposite currents, so both push the same way: a sideways force. Their short sides cross
//! both rows, so their pushes along the track cancel. Like the motor they have the 14.3 mm pitch
//! and three phases commutated with the cart's position (front and rear magnets have opposite
//! polarity too), and the coils under the cart's front and rear halves are driven separately, so
//! the pair gives both sideways force and yaw moment.

use std::f64::consts::PI;

use crate::em::magnet_loops;
use crate::motor::{
    CORNERS, LOOPS, LY, POLARITY, Segment, TAU, motor_segments, mutual_midpoint, phase_currents,
    wrench6,
};

const GUIDE_Z: f64 = -0.00125; // m, second board layer
const GUIDE_TURNS: f64 = 20.0;
/// (ohm) one guidance coil: 20 turns of AWG 24-equivalent copper (84.2 mohm/m) round 2 (TAU + 50 mm).
pub const GUIDE_R: f64 = 20.0 * 2.0 * (TAU + 2.0 * LY) * 0.0842;
/// Coils whose centre lies within this distance ahead of (behind) the cart centre form the front
/// (rear) block.
const BLOCK: f64 = 0.08;

#[derive(Clone, Copy)]
struct GuideCoil {
    x0: f64,
    phase: usize,
}

fn guide_coils(x_from: f64, x_to: f64) -> Vec<GuideCoil> {
    let mut v = Vec::new();
    let mut n = ((x_from - 2.0 * TAU) / (2.0 * TAU)).floor() as i64;
    loop {
        let base = n as f64 * 2.0 * TAU;
        if base > x_to + 2.0 * TAU {
            break;
        }
        for phase in 0..3 {
            let x0 = base + phase as f64 * 2.0 * TAU / 3.0;
            if x0 + TAU >= x_from && x0 <= x_to {
                v.push(GuideCoil { x0, phase });
            }
        }
        n += 1;
    }
    v
}

fn rect(c: &GuideCoil, z: f64) -> [(f64, f64, f64); 4] {
    [
        (c.x0, -LY, z),
        (c.x0 + TAU, -LY, z),
        (c.x0 + TAU, LY, z),
        (c.x0, LY, z),
    ]
}

/// Segments of the guidance coils round a cart at cart_x (and their images), each coil carrying
/// current(coil centre x, phase).
fn guide_segments(cart_x: f64, current: impl Fn(f64, usize) -> f64) -> Vec<Segment> {
    let mut v = Vec::new();
    for c in guide_coils(cart_x - 0.35, cart_x + 0.35) {
        let cur = GUIDE_TURNS * current(c.x0 + 0.5 * TAU, c.phase);
        if cur == 0.0 {
            continue;
        }
        for z in [GUIDE_Z, -GUIDE_Z] {
            let p = rect(&c, z);
            for k in 0..4 {
                let (a, b) = (p[k], p[(k + 1) % 4]);
                v.push((a, (b.0 - a.0, b.1 - a.1, b.2 - a.2), cur));
            }
        }
    }
    v
}

/// Flux linkage of phase `phase`'s guidance coils with the magnets and their images, for a cart at
/// (cart_x, cart_y): coil sides in short pieces, the vector potential at each true midpoint.
fn linkage(cart_x: f64, cart_y: f64, gap: f64, phase: usize, piece: f64) -> f64 {
    let coils: Vec<GuideCoil> = guide_coils(cart_x - 0.35, cart_x + 0.35)
        .into_iter()
        .filter(|c| c.phase == phase)
        .collect();
    let mut psi = 0.0;
    for (k, &(cx, cy)) in CORNERS.iter().enumerate() {
        let (loops, il) = magnet_loops(0.006, 0.004, 1.3, LOOPS, gap);
        let il = il * POLARITY[k];
        let (mx, my) = (cart_x + cx, cart_y + cy);
        for c in &coils {
            let p = rect(c, GUIDE_Z);
            for e in 0..4 {
                let (a, b) = (p[e], p[(e + 1) % 4]);
                let n = ((b.0 - a.0).hypot(b.1 - a.1) / piece).ceil() as usize;
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
                        for zz in [z, -z] {
                            psi += GUIDE_TURNS * il * mutual_midpoint((r, zz, 1.0), p0, p1);
                        }
                    }
                }
            }
        }
    }
    psi
}

/// Guidance wrench tables over one electrical period: per amp of front-block and of rear-block
/// current, [fx, fy, fz, mx, my, mz] about the centre of mass.
pub struct GuideTable {
    pub xs: Vec<f64>,
    pub front: Vec<[f64; 6]>,
    pub rear: Vec<[f64; 6]>,
    pub delta: f64,
    /// (W per A^2) mean copper loss of one block at 1 A.
    pub copper: f64,
}

/// Block current pattern: phase k at electrical angle pi x / tau + delta - 2 pi k / 3, for coils
/// whose centre is in the block.
fn block(cart_x: f64, delta: f64, front: bool) -> impl Fn(f64, usize) -> f64 {
    let i = phase_currents(1.0, cart_x, delta);
    move |xc: f64, phase: usize| {
        let d = xc - cart_x;
        let inside = if front {
            (0.0..BLOCK).contains(&d)
        } else {
            (-BLOCK..0.0).contains(&d)
        };
        if inside { i[phase] } else { 0.0 }
    }
}

impl GuideTable {
    pub fn build(gap: f64, n: usize) -> Self {
        // Phase A alone (all its coils near the cart): F_y(x) ~ R cos(pi x / tau - phi); commutate.
        let (mut s, mut c) = (0.0, 0.0);
        for k in 0..24 {
            let x = 2.0 * TAU * k as f64 / 24.0;
            let segs = guide_segments(x, |_, p| if p == 0 { 1.0 } else { 0.0 });
            let fy = wrench6(x, 0.0, gap, &segs)[1];
            s += fy * (PI * x / TAU).sin();
            c += fy * (PI * x / TAU).cos();
        }
        let delta = -s.atan2(c);
        let xs: Vec<f64> = (0..n).map(|k| 2.0 * TAU * k as f64 / n as f64).collect();
        let front = xs
            .iter()
            .map(|&x| wrench6(x, 0.0, gap, &guide_segments(x, block(x, delta, true))))
            .collect();
        let rear = xs
            .iter()
            .map(|&x| wrench6(x, 0.0, gap, &guide_segments(x, block(x, delta, false))))
            .collect();
        let copper = xs
            .iter()
            .map(|&x| {
                let f = block(x, delta, true);
                guide_coils(x - 0.35, x + 0.35)
                    .iter()
                    .map(|c| f(c.x0 + 0.5 * TAU, c.phase).powi(2) * GUIDE_R)
                    .sum::<f64>()
            })
            .sum::<f64>()
            / n as f64;
        GuideTable {
            xs,
            front,
            rear,
            delta,
            copper,
        }
    }
    /// Per-amp wrenches of the front and rear blocks at cart position x (periodic, linear).
    pub fn at(&self, x: f64) -> ([f64; 6], [f64; 6]) {
        let p = 2.0 * TAU;
        let u = x.rem_euclid(p) / p * self.xs.len() as f64;
        let i = u.floor() as usize % self.xs.len();
        let j = (i + 1) % self.xs.len();
        let t = u - u.floor();
        let mix = |v: &Vec<[f64; 6]>| [0, 1, 2, 3, 4, 5].map(|c| v[i][c] * (1.0 - t) + v[j][c] * t);
        (mix(&self.front), mix(&self.rear))
    }
}

pub fn run() -> GuideTable {
    let gap = 0.005964;
    println!(
        "Guidance coils: {:.1} x {:.0} mm (long sides over the magnet rows), {GUIDE_TURNS:.0} turns ({GUIDE_R:.2} ohm), \
         3 phases at the motor's {:.1} mm pitch, {:.2} mm under the steel; front and rear blocks {:.0} mm long",
        TAU * 1e3,
        2e3 * LY,
        TAU * 1e3,
        -GUIDE_Z * 1e3,
        BLOCK * 1e3
    );
    println!(
        "1. Sideways force from phase A at 1 A: Lorentz on the magnet loops vs coenergy I dPsi/dy (true midpoint)"
    );
    for x in [0.002, 0.006, 0.011] {
        let f = wrench6(
            x,
            0.0,
            gap,
            &guide_segments(x, |_, p| if p == 0 { 1.0 } else { 0.0 }),
        )[1];
        let h = 1e-5;
        for piece in [1e-3, 2.5e-4] {
            let e = (linkage(x, h, gap, 0, piece) - linkage(x, -h, gap, 0, piece)) / (2.0 * h);
            println!(
                "   x {:>4.1} mm, pieces {:.2} mm: Lorentz {f:+.5e} N, coenergy {e:+.5e} N ({:+.1e})",
                x * 1e3,
                piece * 1e3,
                e / f - 1.0
            );
        }
    }
    let t = GuideTable::build(gap, 48);
    let span = |v: Vec<f64>| {
        v.iter()
            .fold((f64::MAX, f64::MIN), |a, &x| (a.0.min(x), a.1.max(x)))
    };
    println!(
        "2. Commutated blocks (offset {:.3} rad), per amp, over one electrical period:",
        t.delta
    );
    for (name, w) in [("front", &t.front), ("rear", &t.rear)] {
        let (fy0, fy1) = span(w.iter().map(|v| v[1]).collect());
        let (mz0, mz1) = span(w.iter().map(|v| v[5]).collect());
        let (fx0, fx1) = span(w.iter().map(|v| v[0]).collect());
        let (fz0, fz1) = span(w.iter().map(|v| v[2]).collect());
        println!(
            "   {name:<5} block: sideways {fy0:+.4} .. {fy1:+.4} N/A, yaw {mz0:+.2e} .. {mz1:+.2e} N m/A; \
             side effects: thrust {fx0:+.1e} .. {fx1:+.1e} N/A, lift {fz0:+.1e} .. {fz1:+.1e} N/A",
        );
    }
    println!("   copper loss per block at 1 A: {:.2} W", t.copper);
    println!("3. Off-centre side effects (cart shifted sideways by 1 mm), per amp:");
    let motor_delta = crate::motor::MotorTable::build(gap, 12).delta;
    for x in [0.0, 0.005, 0.01] {
        let m0 = wrench6(
            x,
            0.0,
            gap,
            &motor_segments(x, phase_currents(1.0, x, motor_delta)),
        );
        let m1 = wrench6(
            x,
            0.001,
            gap,
            &motor_segments(x, phase_currents(1.0, x, motor_delta)),
        );
        let g0 = wrench6(x, 0.0, gap, &guide_segments(x, block(x, t.delta, true)));
        let g1 = wrench6(x, 0.001, gap, &guide_segments(x, block(x, t.delta, true)));
        println!(
            "   x {:>4.1} mm: motor sideways {:+.2e} -> {:+.2e} N/A, thrust {:+.4} -> {:+.4}; guidance front sideways {:+.4} -> {:+.4}",
            x * 1e3,
            m0[1],
            m1[1],
            m0[0],
            m1[0],
            g0[1],
            g1[1]
        );
    }
    t
}
