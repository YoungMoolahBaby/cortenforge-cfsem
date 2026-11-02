//! `link`: inductive power from the track to the maglev cart.
//!
//! A one-turn track loop of Litz wire runs the length of the track in a tray under the cart: its
//! two long conductors (go and return) sit `w` apart. A flat spiral pickup of Litz under the cart's
//! hull spans them. The track is driven at constant current I1 at 85 kHz (series-compensated, so
//! the inverter sees a resistance); the pickup is parallel-tuned (a current source whose voltage
//! the tuning raises by Q2, the usual choice for a pickup feeding a DC bus) into a resistive load
//! R_L standing for the rectifier and charger.
//!
//! Positions follow the maglev cart: rail (steel) underside at z = 0, cart body from -gap to
//! -gap - 30 mm. The steel rail is bounded by images: an ideal iron plane (same current) and a
//! perfect conductor (opposite current), on the rail's underside.

use std::f64::consts::PI;

use cfsem::MU_0;
use nalgebra::Complex;

use crate::loops::{Loop, mutual, self_inductance};

pub const F: f64 = 85e3;
pub const W: f64 = 2.0 * PI * F;
const RHO_CU: f64 = 1.72e-8;
/// Assumed AC / DC resistance of fine-stranded Litz at 85 kHz (to replace with a datasheet value).
const LITZ_AC: f64 = 1.2;

/// Hover gap (m); the cart's other gaps are in `run`.
pub const GAP0: f64 = 0.005964;
/// The pickup winding's plane under the rail: 2 mm under the 30 mm body.
pub fn z_pickup(gap: f64) -> f64 {
    -gap - 0.032
}
/// The track loop runs from x = -0.15 to 1.35 m (the cart travels 0 to 1.2 m).
pub const TRACK_X: (f64, f64) = (-0.15, 1.35);
/// Chosen design (from the sweeps in `run`): conductors 60 mm apart, 8 mm under the pickup at
/// the hover gap; a 130 x 70 mm pickup; loaded Q2 of 8; 5 W into the load.
pub const W_TRACK: f64 = 0.06;
pub const Z_TRACK: f64 = -0.046;
pub const PICKUP: (f64, f64) = (0.13, 0.07);
pub const Q2: f64 = 8.0;
pub const P_LOAD: f64 = 5.0;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Rail {
    Free,
    Iron,
    Conductor,
}
impl Rail {
    fn sign(self) -> f64 {
        match self {
            Rail::Free => 0.0,
            Rail::Iron => 1.0,
            Rail::Conductor => -1.0,
        }
    }
}

/// A coil as a set of series turns (each a one-turn loop), with its round-wire radius and
/// resistance.
#[derive(Clone)]
pub struct Coil {
    pub turns: Vec<Loop>,
    pub wire: f64,
    pub ohm: f64,
}

/// Track: one turn of Litz (2.5 mm^2 copper, bundle radius 1.2 mm), conductors `w` apart in the
/// plane z, from x0 to x1.
pub fn track_span(w: f64, z: f64, x0: f64, x1: f64) -> Coil {
    let t = Loop::rect_by(0.5 * (x0 + x1), 0.0, z, x1 - x0, w, 0.01, 1.0);
    let ohm = t.length() * RHO_CU / 2.5e-6 * LITZ_AC;
    Coil {
        turns: vec![t],
        wire: 1.2e-3,
        ohm,
    }
}

pub fn track(w: f64, z: f64) -> Coil {
    track_span(w, z, TRACK_X.0, TRACK_X.1)
}

/// Pickup: a flat rectangular spiral of `n` turns of Litz (0.5 mm^2 copper, bundle radius
/// 0.55 mm, pitch 1.2 mm), outer size lx x ly, centred under the cart at (x, y).
pub fn pickup(x: f64, y: f64, gap: f64, lx: f64, ly: f64, n: usize, piece: f64) -> Coil {
    let turns: Vec<Loop> = (0..n)
        .map(|k| {
            let s = 2.0 * 1.2e-3 * k as f64;
            Loop::rect_by(x, y, z_pickup(gap), lx - s, ly - s, piece, 1.0)
        })
        .collect();
    let ohm = turns.iter().map(Loop::length).sum::<f64>() * RHO_CU / 0.5e-6 * LITZ_AC;
    Coil {
        turns,
        wire: 0.55e-3,
        ohm,
    }
}

pub fn image(c: &Coil) -> Vec<Loop> {
    c.turns.iter().map(Loop::image).collect()
}

/// (H) mutual inductance, sources `a` (exact segments), targets `b` (3-point Gauss on each piece).
pub fn coil_mutual(a: &[Loop], b: &[Loop]) -> f64 {
    a.iter()
        .map(|s| b.iter().map(|t| mutual(s, t)).sum::<f64>())
        .sum()
}

/// (H) self-inductance: turn self terms (corrected for cfsem's internal term), turn-to-turn
/// mutuals, and the rail image.
pub fn l_self(c: &Coil, rail: Rail) -> f64 {
    let n = c.turns.len();
    let mut l = 0.0;
    for i in 0..n {
        l += self_inductance(&c.turns[i], c.wire);
        for j in 0..n {
            if i != j {
                l += mutual(&c.turns[i], &c.turns[j]);
            }
        }
    }
    if rail != Rail::Free {
        l += rail.sign() * coil_mutual(&image(c), &c.turns);
    }
    l
}

/// (H) track-to-pickup mutual inductance, with the rail image.
pub fn m_link(tr: &Coil, pk: &Coil, rail: Rail) -> f64 {
    let mut m = coil_mutual(&tr.turns, &pk.turns);
    if rail != Rail::Free {
        m += rail.sign() * coil_mutual(&image(tr), &pk.turns);
    }
    m
}

/// The parallel-tuned pickup with its capacitor c2 and load rl fixed, at track current i1 (A rms,
/// phase 0): induced EMF -j w M i1 in series with the coil (L2, R2), then C2 across R_L.
#[derive(Clone, Copy)]
pub struct Operating {
    pub p_load: f64,
    pub p_pickup: f64,
    pub p_track: f64,
    /// (W) eddy loss in the rail: what the track supplies beyond the copper losses and the load.
    pub p_rail: f64,
    pub v_load: f64,
    pub i_coil: f64,
    /// (A rms, complex) the pickup coil's current, phase relative to the track current.
    pub i2: Complex<f64>,
    pub eta: f64,
}

/// As `operate_lossy` with real (lossless) inductances.
pub fn operate(m: f64, l2: f64, r1: f64, r2: f64, c2: f64, rl: f64, i1: f64) -> Operating {
    let c = |v: f64| Complex::new(v, 0.0);
    operate_lossy(c(m), c(l2), c(0.0), r1, r2, c2, rl, i1)
}

/// The link with complex inductances (a lossy rail's image: an imaginary part is a loss): mutual
/// m, pickup l2, and the rail's change to the track's own inductance dl1.
#[allow(clippy::too_many_arguments)]
pub fn operate_lossy(
    m: Complex<f64>,
    l2: Complex<f64>,
    dl1: Complex<f64>,
    r1: f64,
    r2: f64,
    c2: f64,
    rl: f64,
    i1: f64,
) -> Operating {
    let j = Complex::new(0.0, 1.0);
    let zp = Complex::new(rl, 0.0) / (Complex::new(1.0, 0.0) + j * W * c2 * rl);
    let i1c = Complex::new(i1, 0.0);
    let i2 = -j * W * m * i1c / (r2 + j * W * l2 + zp);
    let v = i2 * zp;
    // The track source's voltage, without its copper drop and series capacitor (reactive).
    let v1 = j * W * dl1 * i1c + j * W * m * i2;
    let p_in = (v1 * i1c.conj()).re + i1 * i1 * r1;
    let (p_load, p_pickup, p_track) = (v.norm_sqr() / rl, i2.norm_sqr() * r2, i1 * i1 * r1);
    let p_rail = p_in - p_load - p_pickup - p_track;
    Operating {
        p_load,
        p_pickup,
        p_track,
        p_rail,
        v_load: v.norm(),
        i_coil: i2.norm(),
        i2,
        eta: p_load / p_in,
    }
}

/// The tuned design for an n-turn pickup at the hover gap: (m, l2, track, pickup, c2, rl, i1).
pub struct Design {
    pub m: f64,
    pub l1: f64,
    pub l2: f64,
    pub tr: Coil,
    pub pk: Coil,
    pub c2: f64,
    pub rl: f64,
    pub i1: f64,
    pub op: Operating,
}

pub fn design(n: usize) -> Design {
    let tr = track(W_TRACK, Z_TRACK);
    let pk = pickup(0.6, 0.0, GAP0, PICKUP.0, PICKUP.1, n, 0.005);
    let (m, l1, l2) = (
        m_link(&tr, &pk, Rail::Free),
        l_self(&tr, Rail::Free),
        l_self(&pk, Rail::Free),
    );
    let c2 = 1.0 / (W * W * l2);
    let rl = Q2 * W * l2;
    let p1 = operate(m, l2, tr.ohm, pk.ohm, c2, rl, 1.0).p_load;
    let i1 = (P_LOAD / p1).sqrt();
    let op = operate(m, l2, tr.ohm, pk.ohm, c2, rl, i1);
    Design {
        m,
        l1,
        l2,
        tr,
        pk,
        c2,
        rl,
        i1,
        op,
    }
}

pub fn run() {
    println!(
        "Inductive power for the maglev cart: {:.0} kHz, track loop x {:.2} .. {:.2} m, pickup {:.0} mm under the rail at the hover gap",
        F / 1e3,
        TRACK_X.0,
        TRACK_X.1,
        -z_pickup(GAP0) * 1e3,
    );
    println!(
        "   assumptions: Litz AC/DC resistance {LITZ_AC}; track 2.5 mm^2 copper, pickup 0.5 mm^2 at 1.2 mm pitch; resistive load"
    );

    // 1. One pickup turn over a long track: the flux is the 2D two-wire flux per metre times the
    // turn's length (a long straight pair's field does not vary along it).
    println!(
        "1. One 120 x 60 mm pickup turn at the middle of a 20 m track loop vs the 2D two-wire flux x length"
    );
    let h0 = z_pickup(GAP0) - Z_TRACK;
    for w in [0.04, 0.06, 0.08] {
        let tr = track_span(w, Z_TRACK, -10.0, 10.0);
        let pk = pickup(0.0, 0.0, GAP0, 0.12, 0.06, 1, 0.002);
        let got = m_link(&tr, &pk, Rail::Free);
        // Each wire gives (mu0 / 2 pi) ln(r_far / r_near) through the strip |y| < ly / 2.
        let (ly, hw) = (0.06_f64, 0.5 * w);
        let want = MU_0 / PI * ((0.5 * ly + hw).hypot(h0) / (0.5 * ly - hw).hypot(h0)).ln() * 0.12;
        println!(
            "   w {:.0} mm: cfsem {got:.5e} H, 2D {want:.5e} H ({:+.1e})",
            w * 1e3,
            got / want - 1.0
        );
    }

    // 2. Pickup discretization.
    let tr0 = track(W_TRACK, Z_TRACK);
    for piece in [0.005, 0.0025] {
        let pk = pickup(0.6, 0.0, GAP0, PICKUP.0, PICKUP.1, 12, piece);
        println!(
            "2. 12-turn 130 x 70 mm pickup, pieces {:.1} mm: L2 {:.4} uH, M {:.4} uH",
            piece * 1e3,
            l_self(&pk, Rail::Free) * 1e6,
            m_link(&tr0, &pk, Rail::Free) * 1e6
        );
    }

    // 3. Spacing and height: k^2 Q1 = w M^2 / (L2 R1) sets the track loss (1 / (k^2 Q1 Q2) of the
    // load power) and does not depend on either coil's turns.
    let fom = |w: f64, z: f64| {
        let tr = track(w, z);
        let pk = pickup(0.6, 0.0, GAP0, PICKUP.0, PICKUP.1, 12, 0.005);
        let (m, l1, l2) = (
            m_link(&tr, &pk, Rail::Free),
            l_self(&tr, Rail::Free),
            l_self(&pk, Rail::Free),
        );
        let k2q1 = W * m * m / (l2 * tr.ohm);
        let ql2 = W * l2 / pk.ohm;
        (m, m / (l1 * l2).sqrt(), l1, tr.ohm, k2q1, ql2)
    };
    println!(
        "3. Track conductor spacing (8 mm under the pickup), 12-turn pickup: coupling and the best coil-to-coil efficiency"
    );
    for w in [0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09] {
        let (m, k, l1, r1, k2q1, ql2) = fom(w, Z_TRACK);
        let q = (ql2 / k2q1).sqrt();
        println!(
            "   w {:>2.0} mm: M {:.3} uH, k {:.4}, L1 {:.2} uH, R1 {:.1} mohm, k^2 Q1 {:.3}, pickup Q {:.0}; best at Q2 {:.1}: {:.1}%",
            w * 1e3,
            m * 1e6,
            k,
            l1 * 1e6,
            r1 * 1e3,
            k2q1,
            ql2,
            q,
            100.0 / (1.0 + 1.0 / (k2q1 * q) + q / ql2)
        );
    }
    println!(
        "   track 60 mm wide; height under the pickup (clearance for the tray cover and the cart's 1.4 mm rise on landing):"
    );
    for h in [0.005, 0.008, 0.012, 0.016] {
        let (m, k, _, _, k2q1, ql2) = fom(W_TRACK, z_pickup(GAP0) - h);
        let q = (ql2 / k2q1).sqrt();
        println!(
            "   h {:>2.0} mm: M {:.3} uH, k {:.4}, k^2 Q1 {:.3}; best at Q2 {:.1}: {:.1}%",
            h * 1e3,
            m * 1e6,
            k,
            k2q1,
            q,
            100.0 / (1.0 + 1.0 / (k2q1 * q) + q / ql2)
        );
    }

    // 4. Turns: with the loaded Q2 fixed, turns set the output voltage.
    println!("4. Pickup turns, parallel-tuned at Q2 = {Q2}, {P_LOAD} W into the load:");
    for n in [6, 8, 10, 12, 14] {
        let d = design(n);
        println!(
            "   {n:>2} turns: L2 {:>6.2} uH, M {:.3} uH, R2 {:>3.0} mohm, R_L {:>5.1} ohm, load {:>4.1} V rms; track {:.2} A rms; efficiency {:.1}%",
            d.l2 * 1e6,
            d.m * 1e6,
            d.pk.ohm * 1e3,
            d.rl,
            d.op.v_load,
            d.i1,
            d.op.eta * 100.0
        );
    }
    let n = 10;
    let d = design(n);
    println!(
        "5. Chosen: {n} turns. C2 {:.1} nF ({:.1} V rms, coil {:.2} A rms); track {:.2} A rms, series C1 {:.2} uF, \
         {:.1} V rms across the track's inductance ({:.2} uH); losses: track {:.2} W, pickup {:.2} W; {:.1}% coil to coil",
        d.c2 * 1e9,
        d.op.v_load,
        d.op.i_coil,
        d.i1,
        1.0 / (W * W * d.l1) * 1e6,
        W * d.l1 * d.i1,
        d.l1 * 1e6,
        d.op.p_track,
        d.op.p_pickup,
        d.op.eta * 100.0
    );
    println!(
        "   pickup copper {:.2} m, track {:.2} m",
        d.pk.turns.iter().map(Loop::length).sum::<f64>(),
        d.tr.turns[0].length()
    );

    // 6. Tolerances, capacitor and load fixed at the design, track current fixed.
    println!(
        "6. Load power with C2, R_L and the track current fixed (rail images bound the steel):"
    );
    let case = |x: f64, y: f64, gap: f64, rail: Rail| {
        let pk = pickup(x, y, gap, PICKUP.0, PICKUP.1, n, 0.005);
        let (m, l2) = (m_link(&d.tr, &pk, rail), l_self(&pk, rail));
        (m, l2, operate(m, l2, d.tr.ohm, pk.ohm, d.c2, d.rl, d.i1))
    };
    let row = |label: &str, c: (f64, f64, Operating)| {
        println!(
            "   {label:<38} M {:+.2}%, L2 {:+.3}%: {:.2} W, {:.1} V, {:.1}%",
            (c.0 / d.m - 1.0) * 100.0,
            (c.1 / d.l2 - 1.0) * 100.0,
            c.2.p_load,
            c.2.v_load,
            c.2.eta * 100.0
        );
    };
    for x in [0.0, 0.05, 0.6, 1.15, 1.2] {
        row(
            &format!("along the track, x {x:.2} m"),
            case(x, 0.0, GAP0, Rail::Free),
        );
    }
    for y in [0.001, 0.002] {
        row(
            &format!("sideways {:.0} mm (lips at 2 mm)", y * 1e3),
            case(0.6, y, GAP0, Rail::Free),
        );
    }
    for (label, gap) in [
        ("landed on the skids, gap 4.53 mm", 0.00453),
        ("payload on, mean gap 5.55 mm", 0.00555),
        ("at the station shelf, gap 6.5 mm", 0.0065),
    ] {
        row(label, case(0.6, 0.0, gap, Rail::Free));
    }
    row("rail as ideal iron", case(0.6, 0.0, GAP0, Rail::Iron));
    row(
        "rail as a perfect conductor",
        case(0.6, 0.0, GAP0, Rail::Conductor),
    );
}
