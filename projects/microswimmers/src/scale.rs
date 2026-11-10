//! `scale`: the honest numbers for a person. Whether a swimmer like ours could work inside
//! someone comes down to three things: what a changing field may do to the body (it induces an
//! electric field, which at enough strength fires nerves), how fast a swimmer small enough for
//! the blood vessels can go within that, against the blood's flow, and what coils round a torso
//! take.
//!
//! Stokes flow has no length or time of its own, so the model scales exactly: the drag of moving
//! goes as mu L, of turning as mu L^3, of the coupling between them as mu L^2. A magnet of the
//! same material has a moment M V, V as L^3, so the step-out rate (moment times field over the
//! turning drag) goes as M B / mu whatever the size, and a turn moves the swimmer a fixed
//! fraction of its length. Section 2 checks this against the model run at a thousandth of the
//! size.
//!
//! The body: a field turning at f with strength B0 changes at 2 pi f B0, and the electric field
//! it induces grows with the loops of tissue it can drive current round. For a uniform conductor
//! shaped as an ellipsoid with semi-axes s, a change dB/dt along axis k induces
//! E_i = dB/dt s_i^2 r_j / (s_i^2 + s_j^2), E_j = -dB/dt s_j^2 r_i / (s_i^2 + s_j^2)
//! ((i, j, k) cyclic): its curl is -dB/dt and it runs along the surface, so no charge builds up;
//! that is the quasi-static field exactly. A turning field adds two of these a quarter turn
//! apart. ICNIRP's limits are rms values of this internal field for sinusoids; a turning field
//! is held to the same peak ([e_turning]).
//!
//! The coils: the desk's three pairs with every length (the wire's thickness too) 12.5 times, so
//! the jar's place becomes a cylinder 70 cm across that a patient lies in. The field's evenness
//! is the desk's (it doesn't depend on scale); section 4 checks the scaled coils against the
//! desk's and gives their power, stored energy and drive.

use std::f64::consts::PI;
use std::time::Instant;

use cfsem::MU_0;
use coil_designer::coil::RHO_CU;
use nalgebra::Matrix3;
use rayon::prelude::*;

use crate::check::GLYCEROL_MU;
use crate::coils::{self, C_CU, JAR_R, plane};
use crate::helix::V3;
use crate::rft::swim;
use crate::shape;

const G: f64 = 9.81;
/// Run the field at this fraction of step-out (as in `coils` and `steer`).
const KAPPA: f64 = 0.75;
/// (Pa s, kg/m^3) blood: about three times as viscous as water (Guyton and Hall).
const BLOOD_MU: f64 = 3.0e-3;
const BLOOD_RHO: f64 = 1060.0;
/// (J) k T at body temperature.
const KT: f64 = 1.380649e-23 * 310.0;
/// (m) the body's stand-in, semi-axes: an ellipsoid 36 cm side to side (x), 24 cm front to back
/// (y) and 1.7 m head to foot (z).
const BODY: [f64; 3] = [0.18, 0.12, 0.85];
/// (m) room for the patient: a cylinder 70 cm across along z, like a wide-bore MRI scanner's.
const BORE_R: f64 = 0.35;
/// ICNIRP 2010, Table 2: the rms internal electric field allowed (V/m) at f Hz, in all tissues
/// of head and body (here the trunk) and in the brain and retina. Its limits are a stimulation
/// threshold of 4 V/m divided by 10 (the public) and 5 (at work); the brain's are lower below
/// 1 kHz (400 Hz at work), set by the threshold for phosphenes (flickers of light). The third row
/// is the nerve threshold itself, rising above 3 kHz as the limits do, for the brain too: its
/// drives run above 400 Hz, where ICNIRP has the phosphene and nerve thresholds meet.
struct Limit {
    name: &'static str,
    body: fn(f64) -> f64,
    brain: fn(f64) -> f64,
}

fn public_body(f: f64) -> f64 {
    if f <= 3e3 { 0.4 } else { 1.35e-4 * f }
}
fn public_brain(f: f64) -> f64 {
    if f < 10.0 {
        0.1 / f
    } else if f < 25.0 {
        0.01
    } else if f < 1e3 {
        4e-4 * f
    } else {
        public_body(f)
    }
}
fn work_body(f: f64) -> f64 {
    if f <= 3e3 { 0.8 } else { 2.7e-4 * f }
}
fn work_brain(f: f64) -> f64 {
    if f < 10.0 {
        0.5 / f
    } else if f < 25.0 {
        0.05
    } else if f < 400.0 {
        2e-3 * f
    } else {
        work_body(f)
    }
}
fn nerves(f: f64) -> f64 {
    if f <= 3e3 { 4.0 } else { 1.35e-3 * f }
}

const LIMITS: [Limit; 3] = [
    Limit { name: "public", body: public_body, brain: public_brain },
    Limit { name: "at work", body: work_body, brain: work_brain },
    Limit { name: "nerves fire", body: nerves, brain: nerves },
];
/// (m) the brain's largest field per unit dB/dt: ICNIRP 2010 gives 23-33 mV/m per mT at 50 Hz
/// from detailed body models in a uniform field; the top of that range. A turning field is two
/// such fields a quarter turn apart, so its rms-equivalent (as in [e_turning]) is at most this
/// times 2 pi f B0.
const C_BRAIN: f64 = 0.033 / (2.0 * PI * 50.0 * 1e-3);
/// Guyton and Hall: the total cross-section of each kind of vessel (cm^2). The same flow passes
/// each in turn, 33 cm/s through the aorta's 2.5 cm^2.
const VESSELS: [(&str, f64); 7] = [
    ("aorta", 2.5),
    ("small arteries", 20.0),
    ("arterioles", 40.0),
    ("capillaries", 2500.0),
    ("venules", 250.0),
    ("small veins", 80.0),
    ("venae cavae", 8.0),
];
const AORTA_AREA: f64 = 2.5;
const AORTA_V: f64 = 0.33;
/// The weaker magnet: this fraction of the desk magnet's moment for the same swimmer.
const WEAK: f64 = 0.1;

/// The internal electric field per unit dB/dt along `u`, at `r` in the ellipsoid.
fn e_unit(u: &V3, r: &V3) -> V3 {
    let s2 = [BODY[0].powi(2), BODY[1].powi(2), BODY[2].powi(2)];
    let mut e = V3::zeros();
    for k in 0..3 {
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        let d = s2[i] + s2[j];
        e[i] += u[k] * s2[i] * r[j] / d;
        e[j] -= u[k] * s2[j] * r[i] / d;
    }
    e
}

/// Points on the ellipsoid's surface (where the field is largest: it grows linearly inside).
fn body_surface() -> Vec<V3> {
    let mut v = Vec::new();
    for it in 0..=60 {
        let th = PI * it as f64 / 60.0;
        for ip in 0..120 {
            let ph = 2.0 * PI * ip as f64 / 120.0;
            v.push(V3::new(
                BODY[0] * th.sin() * ph.cos(),
                BODY[1] * th.sin() * ph.sin(),
                BODY[2] * th.cos(),
            ));
        }
    }
    v
}

/// The body's largest field per unit omega B0, the field turning in the plane normal to `n`, as
/// the rms of the sinusoid with the same peak: at each point the field traces an ellipse,
/// E(t) = -a sin wt + b cos wt, and its largest value (the half-major axis) over root 2 is the
/// rms of the field along the direction it points most (ICNIRP's limits are rms values of
/// sinusoids, and peaks are not to be averaged away).
fn e_turning(n: &V3, pts: &[V3]) -> f64 {
    let (e1, e2) = plane(n);
    pts.iter()
        .map(|r| {
            let (a, b) = (e_unit(&e1, r), e_unit(&e2, r));
            let (aa, bb, ab) = (a.norm_squared(), b.norm_squared(), a.dot(&b));
            let major = 0.5 * (aa + bb) + (0.25 * (aa - bb).powi(2) + ab * ab).sqrt();
            (0.5 * major).sqrt()
        })
        .fold(0.0, f64::max)
}

/// The fastest drive a limit allows: the rate f and field B = f / (KAPPA k) (k the step-out rate
/// per tesla) with the trunk's largest field c 2 pi f B (c per unit omega B0) within the limit,
/// and the brain's too if the head is in the field. Both sides grow with f (the field as f^2,
/// the limits at most as f), so the allowed rates run from 0 up to the answer: bisect.
fn safe(lim: &Limit, k: f64, c: f64, head: bool) -> (f64, f64) {
    let ok = |f: f64| {
        let w = 2.0 * PI * f * f / (KAPPA * k);
        c * w <= (lim.body)(f) && (!head || C_BRAIN * w <= (lim.brain)(f))
    };
    // Scan up in small steps first (the limits are piecewise), then bisect the last step.
    let (mut lo, mut f) = (0.0, 1.0);
    while f < 1e6 {
        if ok(f) {
            lo = f;
        }
        f *= 1.01;
    }
    let mut hi = lo * 1.01;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if ok(mid) { lo = mid } else { hi = mid }
    }
    (lo, lo / (KAPPA * k))
}

fn length(x: f64) -> String {
    if x < 1e-6 {
        format!("{:.0} nm", x * 1e9)
    } else if x < 1e-3 {
        format!("{:.1} um", x * 1e6)
    } else if x < 1.0 {
        format!("{:.2} mm", x * 1e3)
    } else {
        format!("{x:.2} m")
    }
}

fn speed(v: f64) -> String {
    if v < 1e-6 {
        format!("{:.0} nm/s", v * 1e9)
    } else if v < 1e-3 {
        format!("{:.1} um/s", v * 1e6)
    } else if v < 1.0 {
        format!("{:.2} mm/s", v * 1e3)
    } else {
        format!("{v:.2} m/s")
    }
}

pub fn run() {
    let t_all = Instant::now();
    let sw = shape::chosen();
    let l0 = sw.length();
    let h = sw.helix;
    let model = sw.model(GLYCEROL_MU);
    let free = swim(&model.r);
    let m0 = sw.moment();
    let ms = sw.mass();
    // Shape constants: lengths moved per turn; the step-out rate per tesla in blood for the desk
    // magnet's material (any size); the translational drag at the desk's size in blood.
    let per_turn = free.u.z * 2.0 * PI / l0;
    let k_blood = m0 / (2.0 * PI * free.torque) * GLYCEROL_MU / BLOOD_MU;
    let a_blood: Matrix3<f64> =
        model.r.fixed_view::<3, 3>(0, 0).into_owned() * (BLOOD_MU / GLYCEROL_MU);
    let mob0 = a_blood.try_inverse().expect("drag invertible")[(2, 2)];
    let width = 2.0 * (h.radius + h.wire).max(sw.head_radius);
    println!(
        "scale: the swimmer `shape` chose (R {:.2} mm, pitch {:.2} mm, {:.1} mm long, {:.2} mm across, {:.1} mg, a {:.1} T magnet) at any size, in blood ({} mPa s, {} kg/m^3)",
        h.radius * 1e3,
        h.pitch * 1e3,
        l0 * 1e3,
        width * 1e3,
        ms.mass * 1e6,
        sw.br,
        BLOOD_MU * 1e3,
        BLOOD_RHO
    );
    println!(
        "   model: {} points at 16 round, extrapolated from 12 and 16 ({:.0} s); axis free",
        model.points, model.seconds
    );
    println!(
        "   a turn moves it {:.4} of its length; step-out in blood {:.0} Hz per mT whatever its size (the desk's magnet material); {:.1}x the swimmer's density",
        per_turn,
        k_blood * 1e-3,
        ms.mass / ms.volume / BLOOD_RHO
    );

    // -----------------------------------------------------------------------------------------
    let pts = body_surface();
    println!(
        "1. What the body allows: the electric field a changing field induces in it (a uniform ellipsoid {:.0} x {:.0} x {:.0} cm)",
        200.0 * BODY[0],
        200.0 * BODY[1],
        200.0 * BODY[2]
    );
    // Check the field: its curl (central differences) against -dB/dt, and how far it leaves the
    // surface (E . normal against |E|), for dB/dt along the axes and a slant.
    let (mut curl_miss, mut leak) = (0.0_f64, 0.0_f64);
    let dl = 1e-4;
    for u in [V3::x(), V3::y(), V3::z(), V3::new(1.0, 2.0, 3.0).normalize()] {
        for r in pts.iter().step_by(97) {
            let r = r * 0.7;
            let d = |ax: usize| {
                let mut e = V3::zeros();
                e[ax] = dl;
                (e_unit(&u, &(r + e)) - e_unit(&u, &(r - e))) / (2.0 * dl)
            };
            let (dx, dy, dz) = (d(0), d(1), d(2));
            let curl = V3::new(dy.z - dz.y, dz.x - dx.z, dx.y - dy.x);
            curl_miss = curl_miss.max((curl + u).norm());
        }
        let top = pts.iter().map(|r| e_unit(&u, r).norm()).fold(0.0, f64::max);
        for r in &pts {
            let n = V3::new(r.x / BODY[0].powi(2), r.y / BODY[1].powi(2), r.z / BODY[2].powi(2))
                .normalize();
            leak = leak.max(e_unit(&u, r).dot(&n).abs() / top);
        }
    }
    println!(
        "   check: curl E against -dB/dt within {curl_miss:.1e}, E . normal on the surface within {leak:.1e} of the largest E"
    );
    let per_mt = |u: V3| {
        pts.iter().map(|r| e_unit(&u, r).norm()).fold(0.0, f64::max) * 2.0 * PI * 50.0 * 1e-3
    };
    println!(
        "   a 50 Hz field, per mT: side to side {:.0} mV/m, front to back {:.0}, head to foot {:.0} at most (ICNIRP 2010: 20-60 in the skin from detailed body models)",
        per_mt(V3::x()) * 1e3,
        per_mt(V3::y()) * 1e3,
        per_mt(V3::z()) * 1e3
    );
    // Every plane the field could turn in (normals over a hemisphere).
    let mut normals = Vec::new();
    for it in 0..=30 {
        let th = 0.5 * PI * it as f64 / 30.0;
        for ip in 0..60 {
            let ph = 2.0 * PI * ip as f64 / 60.0;
            normals.push(V3::new(th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()));
            if it == 0 {
                break;
            }
        }
    }
    let per_plane: Vec<(V3, f64)> = normals
        .par_iter()
        .map(|n| (*n, e_turning(n, &pts)))
        .collect();
    let worst = per_plane.iter().copied().fold((V3::z(), 0.0), |a, b| if b.1 > a.1 { b } else { a });
    let best = per_plane
        .iter()
        .copied()
        .fold((V3::z(), f64::INFINITY), |a, b| if b.1 < a.1 { b } else { a });
    let c = worst.1;
    println!(
        "   a turning field, the largest field (as an rms) per unit 2 pi f B0: {:.3} m at worst (heading ({:.2}, {:.2}, {:.2})), {:.3} m at best (heading ({:.2}, {:.2}, {:.2})); along the body {:.3} m, side to side {:.3}, front to back {:.3}",
        c,
        worst.0.x,
        worst.0.y,
        worst.0.z,
        best.1,
        best.0.x,
        best.0.y,
        best.0.z,
        e_turning(&V3::z(), &pts),
        e_turning(&V3::x(), &pts),
        e_turning(&V3::y(), &pts)
    );
    println!(
        "   so the trunk caps f B0 (any heading, up to 3 kHz; above it the limits rise with frequency):"
    );
    for lim in &LIMITS {
        let e0 = (lim.body)(1e3);
        let fb = e0 / (2.0 * PI * c);
        println!(
            "   {:>12} {e0:>4.1} V/m: f B0 under {:.0} Hz mT (1 mT turning at up to {:.0} Hz, 10 mT at {:.0} Hz)",
            lim.name,
            fb * 1e3,
            fb / 1e-3,
            fb / 10e-3
        );
    }
    println!(
        "   with the head in the field too (the brain {:.0} mV/m per mT at 50 Hz, the top of ICNIRP's 23-33): the public's brain limit caps B0 itself at {:.2} mT from 25 Hz to 1 kHz, at work {:.2} mT from 25 to 400 Hz",
        C_BRAIN * 2.0 * PI * 50.0 * 1e-3 * 1e3,
        4e-4 / (2.0 * PI * C_BRAIN) * 1e3,
        2e-3 / (2.0 * PI * C_BRAIN) * 1e3
    );

    // -----------------------------------------------------------------------------------------
    println!("2. The swimmer at any size");
    let s_chk = 1e-3;
    let t0 = Instant::now();
    let big = sw.model_at(GLYCEROL_MU, 6, 8);
    let small = sw.scaled(s_chk).model_at(BLOOD_MU, 6, 8);
    let mut miss = 0.0_f64;
    for i in 0..6 {
        for j in 0..6 {
            let p = (i >= 3) as i32 + (j >= 3) as i32;
            let pred = big.r[(i, j)] * BLOOD_MU / GLYCEROL_MU * s_chk.powi(1 + p);
            let size = (small.r[(i, i)] * small.r[(j, j)]).sqrt();
            miss = miss.max((small.r[(i, j)] - pred).abs() / size);
        }
    }
    println!(
        "   check: the model at a thousandth of the size in blood against the desk's scaled (mu L, mu L^2, mu L^3): within {miss:.1e} ({:.1} s)",
        t0.elapsed().as_secs_f64()
    );
    // The fastest drive at each limit, for the desk's magnet and a weaker one.
    println!(
        "   the fastest drive each limit allows (at {KAPPA} of step-out, the worst heading), and the speed it gives in lengths per second; the trunk alone in the field, and with the head in it too"
    );
    println!(
        "   {:>12} {:>18} {:>8} {:>8} {:>10} {:>10} {:>8} {:>10}",
        "limit", "magnet", "f Hz", "B0 mT", "lengths/s", "head: f", "B0 mT", "lengths/s"
    );
    // (limit, magnet, head in the field, f, B0)
    let mut drives = Vec::new();
    for lim in &LIMITS {
        for (mag, frac) in [("desk", 1.0), ("tenth", WEAK)] {
            let (f, b) = safe(lim, k_blood * frac, c, false);
            let (fh, bh) = safe(lim, k_blood * frac, c, true);
            println!(
                "   {:>12} {:>18} {f:>8.0} {:>8.3} {:>10.1} {fh:>10.0} {:>8.3} {:>10.1}",
                lim.name,
                if frac == 1.0 { "the desk's" } else { "a tenth as strong" },
                b * 1e3,
                per_turn * f,
                bh * 1e3,
                per_turn * fh
            );
            drives.push((lim.name, mag, false, f, b));
            drives.push((lim.name, mag, true, fh, bh));
        }
    }
    let (_, _, _, f_pub, b_pub) = drives[0];
    println!(
        "   at the public limit with the desk's magnet ({f_pub:.0} Hz, {:.2} mT), by size: speed, sinking nose first, thermal drift in 1 s, the field's grip against heat (m B0 / k T), Reynolds number (the larger of moving and turning; over 1 the model doesn't hold), and the gradient that would pull it as fast as it swims",
        b_pub * 1e3
    );
    println!(
        "   {:>10} {:>12} {:>12} {:>12} {:>12} {:>10} {:>12}",
        "length", "speed", "sinking", "drift 1 s", "m B0 / kT", "Re", "pull T/m"
    );
    let sizes = [2e-6, 5e-6, 10e-6, 30e-6, 100e-6, 300e-6, 1e-3, 3e-3, l0];
    for &l in &sizes {
        let s = l / l0;
        let v = per_turn * l * f_pub;
        let mob = mob0 / s;
        let w = (ms.mass - BLOOD_RHO * ms.volume) * s.powi(3) * G;
        let m = m0 * s.powi(3);
        let re = (BLOOD_RHO * v * l / BLOOD_MU)
            .max(BLOOD_RHO * 2.0 * PI * f_pub * (h.radius * s).powi(2) / BLOOD_MU);
        println!(
            "   {:>10} {:>12} {:>12} {:>12} {:>12.1e} {:>10.2e} {:>12.3e}",
            length(l),
            speed(v),
            speed(mob * w),
            length((2.0 * KT * mob).sqrt()),
            m * b_pub / KT,
            re,
            v / (mob * m)
        );
    }

    // -----------------------------------------------------------------------------------------
    println!(
        "3. Against the blood: the mean flow in each kind of vessel (Guyton and Hall's cross-sections, {:.0} cm/s in the aorta), and the length a swimmer needs to match it",
        AORTA_V * 100.0
    );
    println!("   (* = Reynolds number over 1 at that length: the swimmer would be in flow the model doesn't cover)");
    // The public limit (trunk; head in too) and the nerve threshold, each magnet.
    let cols: Vec<_> = drives
        .iter()
        .filter(|d| d.0 == "public" || (d.0 == "nerves fire" && !d.2))
        .copied()
        .collect();
    print!("   {:>15} {:>11}", "vessel", "mean flow");
    for (name, mag, head, _, _) in &cols {
        let label = if *name == "public" { "public" } else { "nerves" };
        let head = if *head { "+head" } else { "" };
        print!(" {:>17}", format!("{label}{head}, {mag}"));
    }
    println!();
    for (vessel, area) in VESSELS {
        let v = AORTA_V * AORTA_AREA / area;
        print!("   {vessel:>15} {:>11}", speed(v));
        for (_, _, _, f, _) in &cols {
            let l = v / (per_turn * f);
            let re = (BLOOD_RHO * v * l / BLOOD_MU)
                .max(BLOOD_RHO * 2.0 * PI * f * (h.radius * l / l0).powi(2) / BLOOD_MU);
            let mark = if re > 1.0 { "*" } else { " " };
            print!(" {:>16}{mark}", length(l));
        }
        println!();
    }

    // -----------------------------------------------------------------------------------------
    let s = BORE_R / JAR_R;
    let desk = coils::pairs();
    let big = coils::pairs_scaled(s);
    println!(
        "4. The coils for a person: the desk's three pairs, every length {s:.1}x, round a {:.0} cm cylinder (the jar's place) the patient lies in along z",
        200.0 * BORE_R
    );
    let probe = [
        V3::zeros(),
        V3::new(5e-3, 3e-3, 2e-3),
        V3::new(20e-3, -10e-3, 15e-3),
        V3::new(-12e-3, 18e-3, -20e-3),
    ];
    let probe_big: Vec<V3> = probe.iter().map(|p| p * s).collect();
    let (mut e_field, mut e_r, mut e_l, mut e_cu) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let mut rows = Vec::new();
    for (d, b) in desk.iter().zip(&big) {
        let (fd, fb) = (d.field(&probe), b.field(&probe_big));
        for (x, y) in fd.iter().zip(&fb) {
            e_field = e_field.max((y * s - x).norm() / x.norm());
        }
        let (r, l, cu) = (b.coil.resistance(), b.coil.inductance(), b.coil.copper_mass());
        e_r = e_r.max((r * s / d.coil.resistance() - 1.0).abs());
        e_l = e_l.max((l / (d.coil.inductance() * s) - 1.0).abs());
        e_cu = e_cu.max((cu / (d.coil.copper_mass() * s.powi(3)) - 1.0).abs());
        let k = b.field(&[V3::zeros()])[0].dot(&b.axis);
        rows.push((r, l, cu, k));
    }
    println!(
        "   check against the desk's coils: field per amp within {e_field:.1e} of 1/{s:.1} at matching points, resistance {e_r:.1e} of 1/{s:.1}, inductance {e_l:.1e} of {s:.1}x, copper {e_cu:.1e} of {s:.1}^3"
    );
    let surf: Vec<Vec<V3>> = big.iter().map(|p| p.surface()).collect();
    let all = surf.iter().flatten();
    let (mut ext, mut room) = (V3::zeros(), f64::INFINITY);
    for p in all {
        ext = ext.sup(&p.abs());
        room = room.min((p.x * p.x + p.y * p.y).sqrt() - BORE_R);
    }
    let turns = big[0].coil.filaments().len() / 2;
    let wire = big[0].coil.packs[0].wire;
    println!(
        "   the set: {:.2} x {:.2} x {:.2} m outside the windings, {:.1} mm to spare round the 70 cm cylinder; each coil {turns} turns of {:.1} mm copper rod",
        2.0 * ext.x,
        2.0 * ext.y,
        2.0 * ext.z,
        room * 1e3,
        wire.bare * 1e3
    );
    println!(
        "   {:>4} {:>8} {:>10} {:>8} {:>8} {:>8} {:>9} {:>14}",
        "pair", "R m", "mohm", "H", "L/R s", "Cu kg", "mT/A", "kA-turns/mT"
    );
    for (p, &(r, l, cu, k)) in big.iter().zip(&rows) {
        let pack = &p.coil.packs[0];
        println!(
            "   {:>4} {:>8.3} {:>10.2} {:>8.4} {:>8.3} {:>8.0} {:>9.4} {:>14.3}",
            p.name,
            pack.r_mean(),
            r * 1e3,
            l,
            l / r,
            cu,
            k * 1e3,
            1e-3 / k * turns as f64 * 1e-3
        );
    }
    println!(
        "   copper in all: {:.0} kg",
        rows.iter().map(|r| r.2).sum::<f64>()
    );
    // Operating points: each limit's fastest drive, and the desk's top field at the public rate.
    let fb_pub = public_body(1e3) / (2.0 * PI * c);
    let mut points: Vec<(String, f64, f64)> = drives
        .iter()
        .map(|(name, mag, head, f, b)| {
            let head = if *head { " +head" } else { "" };
            (format!("{name}{head}, {mag}"), *f, *b)
        })
        .collect();
    // Drop the head-in rows that match the trunk's.
    points.dedup_by(|b, a| b.1 == a.1 && b.2 == a.2);
    points.push(("5 mT at the public cap".into(), fb_pub / 5e-3, 5e-3));
    println!(
        "   driving them (the worst of heading along x, y or z): a pair's peak current (its two coils in series), mean power, the worst pair's heating with no cooling, stored energy, reactive power (2 pi f times it), the worst pair's peak voltage, copper's skin depth at f"
    );
    println!(
        "   {:>26} {:>7} {:>7} {:>8} {:>9} {:>10} {:>9} {:>10} {:>9} {:>8}",
        "drive", "B0 mT", "f Hz", "peak A", "mean W", "C/min", "stored J", "kVA", "peak kV", "skin mm"
    );
    for (label, f, b0) in &points {
        let w = 2.0 * PI * f;
        let (mut i_pk, mut watts, mut heat, mut stored, mut volts) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for n in [V3::x(), V3::y(), V3::z()] {
            let (e1, e2) = plane(&n);
            let (mut p_n, mut e_n) = (0.0, 0.0);
            for (kk, &(r, l, cu, k)) in rows.iter().enumerate() {
                let i = b0 * (e1[kk].powi(2) + e2[kk].powi(2)).sqrt() / k;
                i_pk = i_pk.max(i);
                p_n += 0.5 * r * i * i;
                e_n += 0.5 * l * i * i;
                heat = heat.max(0.5 * r * i * i / (cu * C_CU) * 60.0);
                volts = volts.max(i * (r * r + (w * l).powi(2)).sqrt());
            }
            watts = watts.max(p_n);
            stored = stored.max(e_n);
        }
        let skin = (RHO_CU / (PI * f * MU_0)).sqrt();
        println!(
            "   {label:>26} {:>7.3} {f:>7.0} {i_pk:>8.1} {watts:>9.2} {heat:>10.2e} {stored:>9.2} {:>10.1} {:>9.2} {:>8.2}",
            b0 * 1e3,
            w * stored * 1e-3,
            volts * 1e-3,
            skin * 1e3
        );
    }
    println!("   ({:.0} s in all)", t_all.elapsed().as_secs_f64());
}
