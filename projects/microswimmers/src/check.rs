//! `check`: the physics before anything is built on it.
//!
//! 1. The regularized Stokeslet against the Stokeslet and incompressibility.
//! 2. A sphere by regularized Stokeslets against Stokes; picks eps / h.
//! 3. Spheroids against the exact solution (Chwang and Wu 1975), the formula itself first.
//! 4. Chirality: a straight wire doesn't swim, a mirrored helix swims backwards.
//! 5. A helix: the regularized-Stokeslet model converging, then resistive force theory against it.
//! 6. CortenForge's own fluid model (sim-core): a sphere, then a helix.
//! 7. cfsem: a Helmholtz pair, and a rotating field from two pairs in quadrature.

use std::f64::consts::PI;
use std::time::Instant;

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_density_circular_filament_scalar;
use nalgebra::Matrix6;

use crate::helix::{Helix, V3, spheroid_points};
use crate::rft::{self, advance, gray_hancock, lighthill, swim};
use crate::{simcore, stokeslet};

/// Pure glycerol at 20 C: viscosity (Pa s) and density (kg/m^3).
pub const GLYCEROL_MU: f64 = 1.41;
pub const GLYCEROL_RHO: f64 = 1261.0;
/// CG tolerance (relative residual) for every regularized-Stokeslet solve.
pub const TOL: f64 = 1e-10;

struct Tally {
    pass: usize,
    fail: usize,
}

impl Tally {
    fn small(&mut self, what: &str, err: f64, tol: f64) {
        let ok = err.abs() <= tol;
        if ok {
            self.pass += 1
        } else {
            self.fail += 1
        }
        println!(
            "   {} {what}: {err:.2e} (tol {tol:.0e})",
            if ok { "ok  " } else { "FAIL" }
        );
    }
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    }
    fn vec(&mut self) -> V3 {
        V3::new(self.next(), self.next(), self.next()) * 2.0
    }
    fn unit(&mut self) -> V3 {
        self.vec().normalize()
    }
}

/// Exact Stokes drag of a prolate spheroid, semi-axis `a` along its axis and `b` across (a > b):
/// force along the axis, across it, and the torque spinning about it, per unit velocity.
pub fn spheroid_exact(a: f64, b: f64, mu: f64) -> (f64, f64, f64) {
    let e = (1.0 - b * b / (a * a)).sqrt();
    let le = 2.0 * e.atanh();
    let e3 = e * e * e;
    (
        16.0 * PI * mu * a * e3 / ((1.0 + e * e) * le - 2.0 * e),
        32.0 * PI * mu * a * e3 / (2.0 * e + (3.0 * e * e - 1.0) * le),
        32.0 * PI * mu * a * b * b * e3 / (3.0 * (2.0 * e - (1.0 - e * e) * le)),
    )
}

fn diag_mean(r: &Matrix6<f64>, i: &[usize]) -> f64 {
    i.iter().map(|&k| r[(k, k)]).sum::<f64>() / i.len() as f64
}

fn asym(r: &Matrix6<f64>) -> f64 {
    (r - r.transpose()).abs().max() / r.abs().max()
}

/// RS resistance of a helix at `n_around` points round the wire, eps = c h.
pub fn rs_helix(h: &Helix, n_around: usize, c: f64, mu: f64) -> (Matrix6<f64>, usize, usize) {
    let (pts, sp) = h.tube_points(n_around);
    let (r, it) = stokeslet::resistance(&pts, c * sp, mu, &h.centre(), TOL);
    (r, pts.len(), it)
}

/// The model's error falls about as the spacing h (section 3), so from 12 and 16 points round
/// the wire (h in the ratio 4 : 3) the fine-spacing limit is R = 4 R16 - 3 R12.
pub fn extrapolate(r12: &Matrix6<f64>, r16: &Matrix6<f64>) -> Matrix6<f64> {
    r16 * 4.0 - r12 * 3.0
}

/// The model of a helix extrapolated from 12 and 16 round; also returns the 16-round matrix
/// and the points at 16 round.
pub fn rs_helix_fine(h: &Helix, c: f64, mu: f64) -> (Matrix6<f64>, Matrix6<f64>, usize) {
    let (r12, _, _) = rs_helix(h, 12, c, mu);
    let (r16, n, _) = rs_helix(h, 16, c, mu);
    (extrapolate(&r12, &r16), r16, n)
}

pub fn run() -> bool {
    let t0 = Instant::now();
    let mut t = Tally { pass: 0, fail: 0 };
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);

    // ---------------------------------------------------------------------------------------
    println!("1. The regularized Stokeslet (Cortez, Fauci and Medovikov 2005)");
    let (mu, eps) = (1.0, 1.0);
    let (mut w100, mut w1000) = (0.0_f64, 0.0_f64);
    for _ in 0..200 {
        let f = rng.unit();
        let dir = rng.unit();
        for (dist, worst) in [(100.0, &mut w100), (1000.0, &mut w1000)] {
            let x = dir * dist;
            let u = stokeslet::velocity(&x, &V3::zeros(), &f, eps, mu);
            let s = (f * (dist * dist) + x * x.dot(&f)) / (8.0 * PI * mu * dist.powi(3));
            *worst = worst.max((u - s).norm() / s.norm());
        }
    }
    t.small("far field vs the Stokeslet at 100 eps, worst of 200", w100, 1e-3);
    t.small("the same at 1000 eps (falls as (eps/r)^2)", w1000, 1e-5);
    let mut worst = 0.0_f64;
    for _ in 0..200 {
        let x = rng.vec() * 1.5;
        let f = rng.unit();
        let dl = 1e-4;
        let mut div = 0.0;
        for i in 0..3 {
            let mut e = V3::zeros();
            e[i] = dl;
            div += (stokeslet::velocity(&(x + e), &V3::zeros(), &f, eps, mu)[i]
                - stokeslet::velocity(&(x - e), &V3::zeros(), &f, eps, mu)[i])
                / (2.0 * dl);
        }
        let scale = stokeslet::velocity(&x, &V3::zeros(), &f, eps, mu).norm() / eps;
        worst = worst.max(div.abs() / scale);
    }
    t.small(
        "divergence by central differences within 3 eps, relative",
        worst,
        1e-6,
    );
    let f = rng.unit();
    let u0 = stokeslet::velocity(&V3::zeros(), &V3::zeros(), &f, eps, mu);
    t.small(
        "velocity at the blob's centre vs f / (4 pi mu eps)",
        u0.norm() * 4.0 * PI * mu * eps - 1.0,
        1e-14,
    );

    // ---------------------------------------------------------------------------------------
    println!(
        "2. A sphere by regularized Stokeslets vs F = 6 pi mu a U, T = 8 pi mu a^3 Omega (a = 1)"
    );
    println!(
        "   error of F / error of T, %; eps = c h, h the point spacing (0.79 a is a tube's 8 round the wire)"
    );
    let cs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75];
    let hs = [0.785, 0.524, 0.393, 0.2, 0.1];
    print!("   {:>6} {:>6}", "h/a", "N");
    for c in cs {
        print!("  {:>13}", format!("c = {c}"));
    }
    println!();
    let mut score = vec![0.0_f64; cs.len()];
    let mut sphere_err = vec![vec![(0.0, 0.0); cs.len()]; hs.len()];
    for (ih, &h) in hs.iter().enumerate() {
        let pts = spheroid_points(1.0, 1.0, h);
        print!("   {h:>6.3} {:>6}", pts.len());
        for (ic, &c) in cs.iter().enumerate() {
            let (r, _) = stokeslet::resistance(&pts, c * h, 1.0, &V3::zeros(), TOL);
            let ef = diag_mean(&r, &[0, 1, 2]) / (6.0 * PI) - 1.0;
            let et = diag_mean(&r, &[3, 4, 5]) / (8.0 * PI) - 1.0;
            sphere_err[ih][ic] = (ef, et);
            score[ic] = score[ic].max(ef.abs().max(et.abs()));
            print!("  {:>13}", format!("{:+.2}/{:+.2}", 100.0 * ef, 100.0 * et));
        }
        println!();
    }
    let best = (0..cs.len())
        .min_by(|&a, &b| score[a].total_cmp(&score[b]))
        .unwrap();
    let c_best = cs[best];
    println!(
        "   c = {c_best} has the least worst error over all spacings ({:.2}%); used from here on",
        100.0 * score[best]
    );
    let (ef, et) = sphere_err[hs.len() - 1][best];
    t.small("sphere drag at h = 0.1 a, relative", ef, 1e-2);
    t.small("sphere torque at h = 0.1 a, relative", et, 1e-2);

    // ---------------------------------------------------------------------------------------
    println!("3. Prolate spheroids: the exact solution (Chwang and Wu 1975), then the model");
    let (fp, fx, tz) = spheroid_exact(1.0, (1.0f64 - 1e-6).sqrt(), 1.0);
    let worst = [fp / (6.0 * PI), fx / (6.0 * PI), tz / (8.0 * PI)]
        .iter()
        .map(|v| (v - 1.0).abs())
        .fold(0.0, f64::max);
    t.small(
        "exact formulas at e = 1e-3 vs the sphere (F along, across, T)",
        worst,
        1e-5,
    );
    let (a, b) = (1000.0, 1.0);
    let (fp, fx, tz) = spheroid_exact(a, b, 1.0);
    let l = (2.0 * a / b).ln();
    let worst = [
        fp / (4.0 * PI * a / (l - 0.5)),
        fx / (8.0 * PI * a / (l + 0.5)),
        tz / (16.0 / 3.0 * PI * a * b * b),
    ]
    .iter()
    .map(|v| (v - 1.0).abs())
    .fold(0.0, f64::max);
    t.small(
        "exact formulas at aspect 1000 vs slender-body theory's leading terms",
        worst,
        1e-4,
    );
    println!(
        "   model vs exact, %: F along / F across / T about the axis (b = 1, eps = {c_best} h)"
    );
    let (mut worst16, mut worst_x) = (0.0_f64, 0.0_f64);
    for aspect in [2.0, 5.0, 10.0] {
        let (fp, fx, tz) = spheroid_exact(aspect, 1.0, 1.0);
        let errs = |r: &Matrix6<f64>| {
            [
                r[(2, 2)] / fp - 1.0,
                diag_mean(r, &[0, 1]) / fx - 1.0,
                r[(5, 5)] / tz - 1.0,
            ]
        };
        let worst_of = |e: &[f64; 3]| e.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let mut by_n = Vec::new();
        for n_around in [8, 12, 16] {
            let h = 2.0 * PI / n_around as f64;
            let pts = spheroid_points(aspect, 1.0, h);
            let (r, it) = stokeslet::resistance(&pts, c_best * h, 1.0, &V3::zeros(), TOL);
            let e = errs(&r);
            if n_around == 16 {
                worst16 = worst16.max(worst_of(&e));
            }
            println!(
                "   aspect {aspect:>4}, {n_around:>2} round ({:>4} points, {it:>4} CG): {:+.2} / {:+.2} / {:+.2}",
                pts.len(),
                100.0 * e[0],
                100.0 * e[1],
                100.0 * e[2]
            );
            by_n.push(r);
        }
        let e = errs(&extrapolate(&by_n[1], &by_n[2]));
        worst_x = worst_x.max(worst_of(&e));
        println!(
            "   aspect {aspect:>4}, extrapolated from 12 and 16 round:  {:+.2} / {:+.2} / {:+.2}",
            100.0 * e[0],
            100.0 * e[1],
            100.0 * e[2]
        );
    }
    println!(
        "   at 16 round alone the worst is {:.2}% (all low); the extrapolated model is used from here on",
        100.0 * worst16
    );
    t.small(
        "worst spheroid error, extrapolated (aspect 2-10)",
        worst_x,
        3e-2,
    );

    // ---------------------------------------------------------------------------------------
    println!(
        "4. Chirality and symmetry (helix R = 1, pitch 2 pi, wire b = 0.2, 2 turns; 12 round)"
    );
    let right = Helix {
        radius: 1.0,
        pitch: 2.0 * PI,
        turns: 2.0,
        wire: 0.2,
        hand: 1.0,
    };
    let left = Helix { hand: -1.0, ..right };
    let straight = Helix {
        radius: 0.0,
        ..right
    };
    let k = lighthill(&right, 1.0);
    let rr = rft::resistance(&right, k, 1.0, 4000);
    let rl = rft::resistance(&left, k, 1.0, 4000);
    let (ar, al) = (advance(&right, &swim(&rr)), advance(&left, &swim(&rl)));
    t.small("RFT: R - R^T, relative", asym(&rr), 1e-12);
    t.small(
        "RFT: mirrored helix advances backwards (sum / size)",
        (ar + al) / ar,
        1e-12,
    );
    let (sr, _, _) = rs_helix(&right, 12, c_best, 1.0);
    let (sl, _, _) = rs_helix(&left, 12, c_best, 1.0);
    let (ar, al) = (advance(&right, &swim(&sr)), advance(&left, &swim(&sl)));
    t.small("model: R - R^T, relative (CG tolerance)", asym(&sr), 1e-6);
    t.small(
        "model: mirrored helix advances backwards (sum / size)",
        (ar + al) / ar,
        1e-6,
    );
    let (ss, _, _) = rs_helix(&straight, 12, c_best, 1.0);
    // Thrust along the axis per unit spin, against the helix's, same wire.
    t.small(
        "model: straight wire spun about its axis, thrust / the helix's",
        ss[(2, 5)] / sr[(2, 5)],
        1e-6,
    );
    println!("   right-handed helix advances {ar:+.4} of its pitch a turn (positive: like a screw)");

    // ---------------------------------------------------------------------------------------
    println!("5. A helix: the model converging, then resistive force theory against it");
    println!(
        "   R = 1, 45 deg, b = 0.2, 2 turns; advance a turn / pitch, torque T / (mu omega R^2 L)"
    );
    let norm_t = |h: &Helix, s: &rft::Swim| s.torque / (h.radius * h.radius * h.length());
    let mut prev = None;
    let mut change = 0.0;
    let mut by_n = Vec::new();
    for n_around in [8, 12, 16] {
        let t1 = Instant::now();
        let (r, n, it) = rs_helix(&right, n_around, c_best, 1.0);
        let s = swim(&r);
        let (adv, tq) = (advance(&right, &s), norm_t(&right, &s));
        println!(
            "   {n_around:>2} round: {n:>5} points, {it:>4} CG, advance {adv:.5}, torque {tq:.5}, wobble {:.1e} ({:.1} s)",
            (s.wobble.0.powi(2) + s.wobble.1.powi(2)).sqrt(),
            t1.elapsed().as_secs_f64()
        );
        if let Some((pa, pt)) = prev {
            let (pa, pt): (f64, f64) = (pa, pt);
            change = f64::max((adv / pa - 1.0).abs(), (tq / pt - 1.0).abs());
        }
        prev = Some((adv, tq));
        by_n.push(r);
    }
    let s = swim(&extrapolate(&by_n[1], &by_n[2]));
    println!(
        "   extrapolated: advance {:.5}, torque {:.5}",
        advance(&right, &s),
        norm_t(&right, &s)
    );
    t.small("model change from 12 to 16 round, relative", change, 2e-2);
    println!(
        "   resistive force theory / the model (extrapolated), 2 turns: Gray-Hancock (GH), Lighthill (L)"
    );
    println!(
        "   {:>5} {:>5} {:>6} {:>9} {:>6} {:>6} {:>9} {:>6} {:>6}",
        "angle", "b/R", "points", "advance", "GH", "L", "torque", "GH", "L"
    );
    for angle in [30.0_f64, 45.0, 60.0] {
        for b in [0.1, 0.2] {
            let h = Helix {
                radius: 1.0,
                pitch: 2.0 * PI / angle.to_radians().tan(),
                turns: 2.0,
                wire: b,
                hand: 1.0,
            };
            let (r, _, n) = rs_helix_fine(&h, c_best, 1.0);
            let s = swim(&r);
            let gh = swim(&rft::resistance(&h, gray_hancock(&h, 1.0), 1.0, 4000));
            let lh = swim(&rft::resistance(&h, lighthill(&h, 1.0), 1.0, 4000));
            let (a0, tq0) = (advance(&h, &s), norm_t(&h, &s));
            println!(
                "   {:>5.1} {b:>5.2} {n:>6} {a0:>9.5} {:>6.3} {:>6.3} {tq0:>9.4} {:>6.3} {:>6.3}",
                h.pitch_angle_deg(),
                advance(&h, &gh) / a0,
                advance(&h, &lh) / a0,
                norm_t(&h, &gh) / tq0,
                norm_t(&h, &lh) / tq0
            );
        }
    }

    // ---------------------------------------------------------------------------------------
    println!("6. CortenForge's fluid model (sim-core 0.9), glycerol ({GLYCEROL_MU} Pa s)");
    let a = 1e-3;
    let (fx, tz) = (1e-4, 1e-8);
    for shape in ["ellipsoid", "none"] {
        let (push, twist) = simcore::sphere(a, GLYCEROL_MU, shape, 1e-6, fx, tz);
        let drag = fx / push.v.x / (6.0 * PI * GLYCEROL_MU * a) - 1.0;
        let rot = tz / twist.w.z / (8.0 * PI * GLYCEROL_MU * a.powi(3)) - 1.0;
        let model = if shape == "ellipsoid" {
            "per-geom ellipsoid model"
        } else {
            "body inertia-box model"
        };
        println!(
            "   sphere, {model}: drag {:+.3}%, torque {:+.3}% against Stokes ({} / {} steps)",
            100.0 * drag,
            100.0 * rot,
            push.steps,
            twist.steps
        );
        if shape == "ellipsoid" {
            t.small("ellipsoid model: sphere drag vs 6 pi mu a", drag, 1e-6);
            t.small("ellipsoid model: sphere torque vs 8 pi mu a^3", rot, 1e-6);
            t.small(
                "xfrc_applied [torque, force]: the pushed sphere doesn't turn, |w| a / |v|",
                push.w.norm() * a / push.v.norm(),
                1e-9,
            );
        }
    }
    let desk = Helix {
        radius: 1.5e-3,
        pitch: 6e-3,
        turns: 2.0,
        wire: 0.3e-3,
        hand: 1.0,
    };
    let (r, _, _) = rs_helix_fine(&desk, c_best, GLYCEROL_MU);
    let s = swim(&r);
    let adv_rs = advance(&desk, &s);
    let omega = 2.0 * PI;
    let torque = s.torque * omega;
    println!(
        "   helix R 1.5 mm, pitch 6 mm, wire 0.3 mm, 2 turns: the model advances {adv_rs:.4} of the pitch a turn"
    );
    for (shape, rho) in [
        ("ellipsoid", 0.0),
        ("none", 0.0),
        ("ellipsoid", GLYCEROL_RHO),
    ] {
        let st = simcore::helix(&desk, 16, GLYCEROL_MU, rho, shape, 1e-6, torque);
        let adv = st.v.z * 2.0 * PI / (st.w.z * desk.pitch);
        println!(
            "   sim-core, {shape:>9}, fluid density {rho:>4}: advance {adv:+.2e} of the pitch a turn, turning drag {:.2}x the model's ({} steps)",
            omega / st.w.z,
            st.steps
        );
        if rho == 0.0 {
            t.small(
                &format!("sim-core ({shape}) helix advance / the model's"),
                adv / adv_rs,
                1e-3,
            );
        }
    }

    // ---------------------------------------------------------------------------------------
    println!("7. cfsem: Helmholtz pairs and a rotating field");
    // A pair along unit axis n: loops of radius a at +-a/2. Field at p.
    let pair = |n: V3, a: f64, ni: f64, p: V3| -> V3 {
        let z = p.dot(&n);
        let rv = p - n * z;
        let r = rv.norm();
        let mut b = V3::zeros();
        for z0 in [-a / 2.0, a / 2.0] {
            let (br, bz) = flux_density_circular_filament_scalar((a, z0, ni), (r, z));
            b += n * bz;
            if r > 0.0 {
                b += rv / r * br;
            }
        }
        b
    };
    let (a1, ni) = (0.05, 100.0);
    let b0 = pair(V3::z(), a1, ni, V3::zeros()).z;
    t.small(
        "Helmholtz centre field vs (4/5)^1.5 mu0 NI / R (cfsem's elliptic fits: ~1e-8)",
        b0 / (0.8f64.powf(1.5) * MU_0 * ni / a1) - 1.0,
        3e-8,
    );
    // x pair at radius a1, y pair nested at 1.2 a1 with its current raised to match.
    let a2 = 1.2 * a1;
    let (mut ripple, mut phase) = (0.0_f64, 0.0_f64);
    for i in 0..16 {
        let wt = 2.0 * PI * i as f64 / 16.0;
        let bv = pair(V3::x(), a1, ni * wt.cos(), V3::zeros())
            + pair(V3::y(), a2, ni * a2 / a1 * wt.sin(), V3::zeros());
        ripple = ripple.max((bv.norm() / b0 - 1.0).abs());
        let mut d = bv.y.atan2(bv.x) - wt;
        d = (d + PI).rem_euclid(2.0 * PI) - PI;
        phase = phase.max(d.abs());
    }
    t.small(
        "two pairs in quadrature: field at the centre, magnitude ripple",
        ripple,
        1e-12,
    );
    t.small("the same: angle - omega t (rad)", phase, 1e-12);

    println!(
        "{} of {} pass ({:.1} s)",
        t.pass,
        t.pass + t.fail,
        t0.elapsed().as_secs_f64()
    );
    t.fail == 0
}
