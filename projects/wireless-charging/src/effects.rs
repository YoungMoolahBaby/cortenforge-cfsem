//! `effects`: side effects of the track-to-cart link.
//!
//! 1. The steel rail at 85 kHz. cfsem gives the two bounds as images (ideal iron, perfect
//!    conductor). Real steel is both, with a skin depth of tens of um: a 2D model of the long
//!    track over a conducting magnetic half-space gives each image term's complex weight
//!    R(k) = (mu_r k - gamma) / (mu_r k + gamma), gamma = sqrt(k^2 + j w mu0 mu_r sigma), per spatial
//!    frequency k across the track. The ratio of that integral to the iron one (R = 1) scales
//!    cfsem's 3D iron image terms; its imaginary part is eddy loss in the rail.
//! 2. 85 kHz voltage induced in the cart's lift coils and in the stator's motor and guidance
//!    coils.
//! 3. Stray field against the ICNIRP 2010 reference levels (Health Physics 99(6), Tables 3 and 4,
//!    3 kHz-10 MHz, unperturbed rms): 27 uT for the general public, 100 uT occupational.

use std::f64::consts::PI;

use cfsem::MU_0;
use nalgebra::Complex;

use crate::link::{
    Coil, GAP0, P_LOAD, PICKUP, Rail, W, W_TRACK, Z_TRACK, coil_mutual, design, image, l_self,
    m_link, operate_lossy, z_pickup,
};
use crate::loops::{Loop, P3, field};

type C = Complex<f64>;

/// Low-carbon steel: conductivity (S/m, assumed) and a range of relative permeabilities (the
/// rail carries the magnets' bias flux; incremental permeability is uncertain).
const SIGMA: f64 = 6e6;
const MU_R: [f64; 3] = [100.0, 300.0, 1000.0];

/// 2D line currents along the track: (y, depth under the rail, weight).
type Lines = Vec<(f64, f64, f64)>;

fn r_coef(k: f64, mu_r: f64) -> C {
    let g = C::new(k * k, W * MU_0 * mu_r * SIGMA).sqrt();
    (mu_r * k - g) / (mu_r * k + g)
}

/// Reflected 2D linkage of `tgt` from `src` (units mu0 / 2pi per metre per amp): the sum over
/// line pairs of w_s w_t integral R(k) cos(k dy) e^{-k (d_s + d_t)} / k dk. Both sets carry
/// zero net current, so the sum vanishes as k -> 0 and the integral converges there.
fn reflected(src: &Lines, tgt: &Lines, r: &dyn Fn(f64) -> C) -> C {
    let dmin = src.iter().map(|s| s.1).fold(f64::MAX, f64::min)
        + tgt.iter().map(|t| t.1).fold(f64::MAX, f64::min);
    let (kmax, n) = (60.0 / dmin, 60_000);
    let h = kmax / n as f64;
    let mut sum = C::new(0.0, 0.0);
    for i in 0..n {
        let k = (i as f64 + 0.5) * h;
        let mut s = 0.0;
        for &(ys, ds, ws) in src {
            for &(yt, dt, wt) in tgt {
                s += ws * wt * (k * (yt - ys)).cos() * (-k * (ds + dt)).exp();
            }
        }
        sum += r(k) * s / k;
    }
    sum * h
}

/// Closed forms (same units): an image above the plane (R = 1), and the direct linkage.
fn image_closed(src: &Lines, tgt: &Lines) -> f64 {
    src.iter()
        .flat_map(|&(ys, ds, ws)| {
            tgt.iter()
                .map(move |&(yt, dt, wt)| -ws * wt * (yt - ys).hypot(ds + dt).ln())
        })
        .sum()
}
fn direct_closed(src: &Lines, tgt: &Lines) -> f64 {
    src.iter()
        .flat_map(|&(ys, ds, ws)| {
            tgt.iter()
                .map(move |&(yt, dt, wt)| -ws * wt * (yt - ys).hypot(ds - dt).ln())
        })
        .sum()
}

fn pair(half: f64, d: f64) -> Lines {
    vec![(half, d, 1.0), (-half, d, -1.0)]
}

/// Lift coil: 80 turns of AWG 28 (0.381 mm over the enamel), 8 layers from r 6.5 mm, 10 turns
/// axially, top at the gap; as one-turn 48-gons.
fn lift_coil(cx: f64, cy: f64, gap: f64) -> Vec<Loop> {
    let p = 0.321e-3 + 0.06e-3;
    let mut v = Vec::new();
    for l in 0..8 {
        for k in 0..10 {
            v.push(Loop::circle(
                (cx, cy, -gap - (k as f64 + 0.5) * p),
                6.5e-3 + (l as f64 + 0.5) * p,
                48,
                1.0,
            ));
        }
    }
    v
}

/// Flux linkage of `tgt` per amp in `src`, with a rail image weight (1 iron, -1 conductor).
fn linkage(src: &Coil, tgt: &[Loop], rail: f64) -> f64 {
    let mut m = coil_mutual(&src.turns, tgt);
    if rail != 0.0 {
        m += rail * coil_mutual(&image(src), tgt);
    }
    m
}

pub fn run() {
    let d = design(10);
    let tr = &d.tr;
    let pk = &d.pk;
    let (dt, dp) = (-Z_TRACK, -z_pickup(GAP0));
    let track2 = pair(0.5 * W_TRACK, dt);
    let pick2: Lines = (0..10)
        .flat_map(|k| pair(0.5 * (PICKUP.1 - 2.4e-3 * k as f64), dp))
        .collect();

    println!(
        "1. The steel rail at 85 kHz (2D track over a conducting magnetic half-space, sigma {:.0} MS/m)",
        SIGMA / 1e6
    );
    // Check the spectral integral against the closed-form image (R = 1) and its R = -1 mirror.
    for (name, s, t) in [
        ("track to pickup", &track2, &pick2),
        ("track to itself", &track2, &track2),
        ("pickup to itself", &pick2, &pick2),
    ] {
        let got = reflected(s, t, &|_| C::new(1.0, 0.0));
        let want = image_closed(s, t);
        let neg = reflected(s, t, &|_| C::new(-1.0, 0.0));
        println!(
            "   check, {name}: spectral iron image {:.6e} vs closed form {:.6e} ({:+.1e}); conductor {:+.1e} off its negative",
            got.re,
            want,
            got.re / want - 1.0,
            neg.re / want + 1.0
        );
    }
    let skin = |mu_r: f64| (2.0 / (W * MU_0 * mu_r * SIGMA)).sqrt();
    // cfsem's 3D iron image terms.
    let dm_iron = m_link(tr, pk, Rail::Iron) - d.m;
    let dl2_iron = l_self(pk, Rail::Iron) - d.l2;
    let dl1_iron = l_self(tr, Rail::Iron) - d.l1;
    println!(
        "   cfsem 3D iron image terms: M {:+.4} uH, L2 {:+.4} uH, L1 {:+.4} uH; 2D x length for M: {:+.4} uH",
        dm_iron * 1e6,
        dl2_iron * 1e6,
        dl1_iron * 1e6,
        // Each pickup turn: the 2D image linkage per metre times the turn's length along the track.
        MU_0 / (2.0 * PI)
            * (0..10)
                .map(
                    |k| image_closed(&track2, &pair(0.5 * (PICKUP.1 - 2.4e-3 * k as f64), dp))
                        * (PICKUP.0 - 2.4e-3 * k as f64)
                )
                .sum::<f64>()
            * 1e6
    );
    let mut need = 0.0_f64;
    for mu_r in MU_R {
        let rr = |k: f64| r_coef(k, mu_r);
        let rho = |s: &Lines, t: &Lines| reflected(s, t, &rr) / image_closed(s, t);
        let (rm, r2, r1) = (
            rho(&track2, &pick2),
            rho(&pick2, &pick2),
            rho(&track2, &track2),
        );
        let m = C::new(d.m, 0.0) + rm * dm_iron;
        let l2 = C::new(d.l2, 0.0) + r2 * dl2_iron;
        let dl1 = r1 * dl1_iron;
        let op = operate_lossy(m, l2, dl1, tr.ohm, pk.ohm, d.c2, d.rl, d.i1);
        let i1 = d.i1 * (P_LOAD / op.p_load).sqrt();
        let op5 = operate_lossy(m, l2, dl1, tr.ohm, pk.ohm, d.c2, d.rl, i1);
        need = need.max(i1);
        println!(
            "   mu_r {mu_r:>4.0} (skin depth {:.0} um): image weights M {:+.3}{:+.3}j, L2 {:+.3}{:+.3}j, L1 {:+.3}{:+.3}j; \
             M {:+.2}%, L2 {:+.2}%; at {:.2} A: {:.2} W",
            skin(mu_r) * 1e6,
            rm.re,
            rm.im,
            r2.re,
            r2.im,
            r1.re,
            r1.im,
            (m.re / d.m - 1.0) * 100.0,
            (l2.re / d.l2 - 1.0) * 100.0,
            d.i1,
            op.p_load
        );
        println!(
            "        for {P_LOAD} W: track {i1:.2} A; losses track {:.2} W, pickup {:.2} W, rail {:.2} W; {:.1}% track to load",
            op5.p_track,
            op5.p_pickup,
            op5.p_rail,
            op5.eta * 100.0
        );
    }
    println!("   sized for the worst case: track {need:.2} A rms");

    // The pickup current's phasor at that track current (free-space circuit; the rail changes it
    // by a few % and is bounded below).
    let i1 = need;
    let op = operate_lossy(
        C::new(d.m, 0.0),
        C::new(d.l2, 0.0),
        C::new(0.0, 0.0),
        tr.ohm,
        pk.ohm,
        d.c2,
        d.rl,
        i1,
    );
    let i2 = op.i2;
    println!(
        "2. Induced 85 kHz voltage at track {i1:.2} A and pickup {:.2} A rms ({:+.0} deg), rail: free / ideal iron / perfect conductor",
        i2.norm(),
        i2.arg().to_degrees()
    );
    let volts = |tgt: &[Loop], rail: f64| {
        (C::new(linkage(tr, tgt, rail) * i1, 0.0) + linkage(pk, tgt, rail) * i2).norm() * W
    };
    let l_lift = 0.122e-3; // H, near the rail (maglev-train `corner`)
    for (k, &(cx, cy)) in [
        (-0.05, -0.025),
        (0.05, -0.025),
        (0.05, 0.025),
        (-0.05, 0.025),
    ]
    .iter()
    .enumerate()
    {
        let coil = lift_coil(0.6 + cx, cy, GAP0);
        let v = [0.0, 1.0, -1.0].map(|r| volts(&coil, r));
        println!(
            "   lift coil {k} ({:+.0}, {:+.0}) mm: {:.1} / {:.1} / {:.1} mV rms; through its 0.122 mH: {:.2} / {:.2} / {:.2} mA",
            cx * 1e3,
            cy * 1e3,
            v[0] * 1e3,
            v[1] * 1e3,
            v[2] * 1e3,
            v[0] / (W * l_lift) * 1e3,
            v[1] / (W * l_lift) * 1e3,
            v[2] / (W * l_lift) * 1e3
        );
    }
    // Stator coils lie on the rail's face, where the rail decides the normal field: iron doubles
    // it, a conductor cancels it. Realistic factor from the 2D model (track source).
    let tau = 0.1 / 7.0;
    for (name, coil, lines) in [
        (
            "guidance coil under the cart",
            Loop::rect_by(0.6, 0.0, -0.00125, tau, 0.05, 0.002, 20.0),
            pair(0.025, 0.00125),
        ),
        (
            "guidance coil 0.3 m from it",
            Loop::rect_by(0.3, 0.0, -0.00125, tau, 0.05, 0.002, 20.0),
            pair(0.025, 0.00125),
        ),
        (
            "motor coil under the cart",
            Loop::rect_by(0.6, 0.025, -0.00075, tau, 0.024, 0.002, 20.0),
            vec![(0.037, 0.00075, 1.0), (0.013, 0.00075, -1.0)],
        ),
    ] {
        let tgt = [coil];
        let v = [0.0, 1.0, -1.0].map(|r| volts(&tgt, r));
        let inc = direct_closed(&track2, &lines);
        let f: Vec<f64> = MU_R
            .iter()
            .map(|&mu_r| {
                ((C::new(inc, 0.0) + reflected(&track2, &lines, &|k| r_coef(k, mu_r))) / inc).norm()
            })
            .collect();
        println!(
            "   {name:<30} {:.0} / {:.0} / {:.1} mV rms; steel factor on the free value (2D) {:.3} / {:.3} / {:.3} for mu_r 100 / 300 / 1000: {:.1} to {:.1} mV",
            v[0] * 1e3,
            v[1] * 1e3,
            v[2] * 1e3,
            f[0],
            f[1],
            f[2],
            v[0] * f.iter().cloned().fold(f64::MAX, f64::min) * 1e3,
            v[0] * f.iter().cloned().fold(0.0, f64::max) * 1e3
        );
    }

    println!(
        "3. Stray field (rms) at track {i1:.2} A with the cart over the point (x 0.6 m) and away from it (x 0.2 m), worst of the rail bounds"
    );
    let b_at = |pts: &[P3]| -> Vec<f64> {
        let mut worst = vec![0.0_f64; pts.len()];
        for rail in [0.0, 1.0, -1.0] {
            let tr_img = image(tr);
            let pk_img = image(pk);
            let mut srcs_t: Vec<(&Loop, f64)> = tr.turns.iter().map(|l| (l, 1.0)).collect();
            let mut srcs_p: Vec<(&Loop, f64)> = pk.turns.iter().map(|l| (l, 1.0)).collect();
            if rail != 0.0 {
                srcs_t.extend(tr_img.iter().map(|l| (l, rail)));
                srcs_p.extend(pk_img.iter().map(|l| (l, rail)));
            }
            let bt = field(&srcs_t, pts);
            let bp = field(&srcs_p, pts);
            for (k, w) in worst.iter_mut().enumerate() {
                let c = |a: f64, b: f64| C::new(a * i1, 0.0) + i2 * b;
                let b = (c(bt[k].0, bp[k].0).norm_sqr()
                    + c(bt[k].1, bp[k].1).norm_sqr()
                    + c(bt[k].2, bp[k].2).norm_sqr())
                .sqrt();
                *w = w.max(b);
            }
        }
        worst
    };
    let dists: Vec<f64> = (1..=60).map(|k| 0.01 * k as f64).collect();
    for x in [0.6, 0.2] {
        for (dir, pts) in [
            (
                "below the track",
                dists
                    .iter()
                    .map(|&s| (x, 0.0, Z_TRACK - s))
                    .collect::<Vec<P3>>(),
            ),
            (
                "beside it",
                dists
                    .iter()
                    .map(|&s| (x, 0.5 * W_TRACK + s, Z_TRACK))
                    .collect::<Vec<P3>>(),
            ),
        ] {
            let b = b_at(&pts);
            let reach = |lim: f64| {
                dists
                    .iter()
                    .zip(&b)
                    .filter(|(_, b)| **b > lim)
                    .map(|(s, _)| *s)
                    .fold(0.0, f64::max)
            };
            let at = |s: f64| {
                b[dists
                    .iter()
                    .position(|&v| (v - s).abs() < 1e-9)
                    .expect("grid")]
                    * 1e6
            };
            println!(
                "   x {x:.1} m, {dir:<16}: {:.0} uT at 5 cm, {:.1} at 10 cm, {:.1} at 20 cm, {:.2} at 50 cm; above 100 uT to {:.0} cm, above 27 uT to {:.0} cm",
                at(0.05),
                at(0.10),
                at(0.20),
                at(0.50),
                reach(100e-6) * 100.0,
                reach(27e-6) * 100.0
            );
        }
    }
}
