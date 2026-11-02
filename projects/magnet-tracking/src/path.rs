//! `path`: tracking a moving magnet frame by frame.
//!
//! The 4 x 4 MMC5983MA board, calibrated (noise only), the magnet's moment known. The magnet runs
//! once round a 30 mm circle with its axis along its travel (tilted 20 deg up, like a capsule),
//! read at the sensor's four rate settings (MEMSIC MMC5983MA Rev. A: 0.4 / 0.6 / 0.8 / 1.2 mG at
//! 50 / 100 / 225 / 580 Hz). Each frame is solved from the last one's pose.
//!
//! Sensors read at once (each triggered together over SPI) see one pose; read one after another
//! across the frame, each sees the magnet where it was at its own moment. A static solver then
//! mixes poses; a moving one also fits the magnet's velocity and turn rate within the frame from
//! each sensor's sample time.

use std::f64::consts::PI;
use std::time::Instant;

use nalgebra::DVector;

use crate::calib::lm;
use crate::layout::{SENSORS, weights};
use crate::magnet::V3;
use crate::solve::{EARTH, MAG, Rng};
use crate::tracker::{Array, Pose};

/// Rate settings: (output rate Hz, noise per axis T).
const RATES: [(f64, f64); 4] = [
    (50.0, 0.04e-6),
    (100.0, 0.06e-6),
    (225.0, 0.08e-6),
    (580.0, 0.12e-6),
];
const RADIUS: f64 = 0.03;

/// The truth at time t for speed v and height z: position and axis.
fn truth(t: f64, v: f64, z: f64) -> (V3, V3) {
    let a = v * t / RADIUS;
    (
        V3::new(RADIUS * a.cos(), RADIUS * a.sin(), z),
        V3::new(-a.sin(), a.cos(), 20f64.to_radians().tan()).normalize(),
    )
}

fn tangents(u: &V3) -> (V3, V3) {
    let a = if u.x.abs() < 0.9 { V3::x() } else { V3::y() };
    let e1 = u.cross(&a).normalize();
    (e1, u.cross(&e1))
}

/// A frame's readings: sensor k sampled at offset `tau[k]` from the frame time.
fn readings(
    array: &Array,
    t: f64,
    tau: &[f64],
    v: f64,
    z: f64,
    noise: f64,
    rng: &mut Rng,
) -> DVector<f64> {
    let mut out = DVector::zeros(3 * array.sensors.len());
    for (k, s) in array.sensors.iter().enumerate() {
        let (p, u) = truth(t + tau[k], v, z);
        let b = MAG.exact(&p, &u, s) + EARTH;
        for a in 0..3 {
            out[3 * k + a] = b[a] + noise * rng.normal();
        }
    }
    out
}

/// How a frame is solved.
#[derive(Clone, Copy, PartialEq)]
enum Solver {
    /// One pose for every sensor (8 unknowns: position, axis, background).
    Static,
    /// The same 8, each sensor's pose moved by a motion predicted from the last frames.
    Predicted,
    /// Velocity and axis rate fitted too (13 unknowns), each sensor at its own time.
    Moving,
}

/// A frame's weighted residuals for unknowns x around `init`: position (3), axis (2, tangent),
/// background (3), and with `motion` the velocity (3) and axis rate (2, tangent); otherwise the
/// motion is held at (`vel`, `wv`).
#[allow(clippy::too_many_arguments)]
fn frame_resid(
    x: &DVector<f64>,
    meas: &DVector<f64>,
    array: &Array,
    tau: &[f64],
    w: &DVector<f64>,
    init: &Pose,
    motion: bool,
    vel: V3,
    wv: V3,
) -> DVector<f64> {
    let (e1, e2) = tangents(&init.u);
    let p = init.p + V3::new(x[0], x[1], x[2]);
    let u = (init.u + e1 * x[3] + e2 * x[4]).normalize();
    let g = init.g + V3::new(x[5], x[6], x[7]);
    let (v, om) = if motion {
        (V3::new(x[8], x[9], x[10]), e1 * x[11] + e2 * x[12])
    } else {
        (vel, wv)
    };
    let mut r = DVector::zeros(meas.len());
    for (k, s) in array.sensors.iter().enumerate() {
        let pk = p + v * tau[k];
        let uk = (u + om * tau[k]).normalize();
        let b = crate::tracker::dipole(&pk, &(uk * init.m), s) + g;
        for a in 0..3 {
            r[3 * k + a] = (b[a] - meas[3 * k + a]) * w[3 * k + a];
        }
    }
    r
}

/// Solve a frame; returns the pose at the frame time, its velocity and axis rate.
#[allow(clippy::too_many_arguments)]
fn solve_frame(
    meas: &DVector<f64>,
    array: &Array,
    tau: &[f64],
    w: &DVector<f64>,
    init: &Pose,
    motion: bool,
    vel: V3,
    wv: V3,
) -> (Pose, V3, V3) {
    let (e1, e2) = tangents(&init.u);
    let n = if motion { 13 } else { 8 };
    let mut x0 = DVector::zeros(n);
    if motion {
        x0[8] = vel.x;
        x0[9] = vel.y;
        x0[10] = vel.z;
        x0[11] = wv.dot(&e1);
        x0[12] = wv.dot(&e2);
    }
    let h = [
        1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-10, 1e-10, 1e-10, 1e-5, 1e-5, 1e-5, 1e-4, 1e-4,
    ];
    let x = lm(
        &|x: &DVector<f64>| frame_resid(x, meas, array, tau, w, init, motion, vel, wv),
        x0,
        &h[..n],
    );
    let p = init.p + V3::new(x[0], x[1], x[2]);
    let u = (init.u + e1 * x[3] + e2 * x[4]).normalize();
    let g = init.g + V3::new(x[5], x[6], x[7]);
    let (v, om) = if motion {
        (V3::new(x[8], x[9], x[10]), e1 * x[11] + e2 * x[12])
    } else {
        (vel, wv)
    };
    (Pose { p, u, m: init.m, g }, v, om)
}

/// The position bound (root of the covariance trace) of a frame's solver at the true pose and
/// motion.
#[allow(clippy::too_many_arguments)]
fn bound(
    array: &Array,
    tau: &[f64],
    w: &DVector<f64>,
    truth: &Pose,
    motion: bool,
    vel: V3,
    wv: V3,
) -> f64 {
    let n = if motion { 13 } else { 8 };
    let (e1, e2) = tangents(&truth.u);
    let mut x0 = DVector::zeros(n);
    if motion {
        x0[8] = vel.x;
        x0[9] = vel.y;
        x0[10] = vel.z;
        x0[11] = wv.dot(&e1);
        x0[12] = wv.dot(&e2);
    }
    let zero = DVector::zeros(w.len());
    let h = [
        1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-10, 1e-10, 1e-10, 1e-5, 1e-5, 1e-5, 1e-4, 1e-4,
    ];
    let mut j = nalgebra::DMatrix::zeros(w.len(), n);
    for i in 0..n {
        let (mut a, mut b) = (x0.clone(), x0.clone());
        a[i] += h[i];
        b[i] -= h[i];
        let col = (frame_resid(&a, &zero, array, tau, w, truth, motion, vel, wv)
            - frame_resid(&b, &zero, array, tau, w, truth, motion, vel, wv))
            / (2.0 * h[i]);
        j.set_column(i, &col);
    }
    let c = (j.transpose() * j).try_inverse().expect("identifiable");
    (c[(0, 0)] + c[(1, 1)] + c[(2, 2)]).sqrt()
}

struct Run {
    rms: f64,
    axis_rms: f64,
    worst: f64,
    crb_rms: f64,
    over_1mm: usize,
    frames: usize,
    us_per_frame: f64,
}

/// One trip round the circle.
fn trip(
    array: &Array,
    v: f64,
    z: f64,
    rate: f64,
    noise: f64,
    sequential: bool,
    solver: Solver,
) -> Run {
    let sensor = SENSORS[0];
    let n = array.sensors.len();
    let frame = 1.0 / rate;
    let tau: Vec<f64> = (0..n)
        .map(|k| {
            if sequential {
                ((k as f64 + 0.5) / n as f64 - 0.5) * frame
            } else {
                0.0
            }
        })
        .collect();
    let frames = ((2.0 * PI * RADIUS / v) * rate).round() as usize;
    let m0 = MAG.moment();
    let mut rng = Rng(0x1234_5678_9abc_def1 ^ (v * 1e3) as u64 ^ (rate as u64) << 20);
    let (p0, u0) = truth(0.0, v, z);
    let mut est = Pose {
        p: p0 + V3::new(1e-3, -1e-3, 1e-3),
        u: u0,
        m: m0,
        g: V3::zeros(),
    };
    let (mut vel, mut wv) = (V3::zeros(), V3::zeros());
    let (mut se, mut sa, mut sc, mut worst, mut over) = (0.0, 0.0, 0.0, 0.0_f64, 0);
    let mut solve_time = 0.0;
    for f in 0..frames {
        let t = f as f64 * frame;
        let meas = readings(array, t, &tau, v, z, noise, &mut rng);
        let (mut w, _) = weights(&sensor, &meas);
        w.iter_mut().for_each(|x| {
            if *x > 0.0 {
                *x = 1.0 / noise;
            }
        });
        let prev = est;
        let t0 = Instant::now();
        let (pose, nv, nw) = match solver {
            Solver::Static => solve_frame(
                &meas,
                array,
                &tau,
                &w,
                &est,
                false,
                V3::zeros(),
                V3::zeros(),
            ),
            Solver::Predicted => solve_frame(&meas, array, &tau, &w, &est, false, vel, wv),
            Solver::Moving => solve_frame(&meas, array, &tau, &w, &est, true, vel, wv),
        };
        solve_time += t0.elapsed().as_secs_f64();
        est = pose;
        match solver {
            Solver::Static => {}
            // The prediction: an average of the last frames' motion (half weight on the newest).
            Solver::Predicted => {
                vel = vel * 0.5 + (est.p - prev.p) * (0.5 * rate);
                wv = wv * 0.5 + (est.u - prev.u) * (0.5 * rate);
            }
            Solver::Moving => {
                vel = nv;
                wv = nw;
            }
        }
        let (p, u) = truth(t, v, z);
        let e = (est.p - p).norm();
        if f >= 5 {
            let truth_pose = Pose {
                p,
                u,
                m: m0,
                g: EARTH,
            };
            let tv = (truth(t + 1e-6, v, z).0 - truth(t - 1e-6, v, z).0) / 2e-6;
            let tw = (truth(t + 1e-6, v, z).1 - truth(t - 1e-6, v, z).1) / 2e-6;
            let b = bound(
                array,
                &tau,
                &w,
                &truth_pose,
                solver == Solver::Moving,
                tv,
                tw,
            );
            se += e * e;
            sa += u.angle(&est.u).powi(2);
            sc += b * b;
            worst = worst.max(e);
            if e > 1e-3 {
                over += 1;
            }
        }
    }
    let k = (frames - 5) as f64;
    Run {
        rms: (se / k).sqrt(),
        axis_rms: (sa / k).sqrt(),
        worst,
        crb_rms: (sc / k).sqrt(),
        over_1mm: over,
        frames,
        us_per_frame: solve_time * 1e6 / frames as f64,
    }
}

pub fn run() {
    let array = Array::grid(4, 4, 0.03);
    println!(
        "Tracking round a {:.0} mm circle (axis along the travel, 20 deg up), 4 x 4 MMC5983MA at 30 mm, calibrated, moment \
         known; one trip each; error RMS over the trip, against the solver's own Cramer-Rao bound along it",
        RADIUS * 1e3
    );
    let fmt = |x: f64| {
        if x < 1e-3 {
            format!("{:.1} um", x * 1e6)
        } else {
            format!("{:.2} mm", x * 1e3)
        }
    };
    println!("1. Sensors read at once, static solver, 50 mm/s:");
    for z in [0.04, 0.08] {
        for (rate, noise) in RATES {
            let r = trip(&array, 0.05, z, rate, noise, false, Solver::Static);
            println!(
                "   {:.0} mm up, {:>3.0} Hz ({:.2} uT): {} RMS (bound {}), axis {:.3} deg, worst {}; {} frames, {:.0} us a frame",
                z * 1e3,
                rate,
                noise * 1e6,
                fmt(r.rms),
                fmt(r.crb_rms),
                r.axis_rms.to_degrees(),
                fmt(r.worst),
                r.frames,
                r.us_per_frame
            );
        }
    }
    println!(
        "2. Sensors read one after another across the frame, 40 mm up: RMS (bound), worst, frames over 1 mm, solve time"
    );
    for (rate, noise) in [RATES[0], RATES[2]] {
        for v in [0.01, 0.05, 0.2] {
            println!(
                "   {:>3.0} Hz, {:>3.0} mm/s ({:.2} mm a frame):",
                rate,
                v * 1e3,
                v / rate * 1e3
            );
            for (name, solver) in [
                ("static", Solver::Static),
                ("predicted motion", Solver::Predicted),
                ("fitted motion", Solver::Moving),
            ] {
                let r = trip(&array, v, 0.04, rate, noise, true, solver);
                println!(
                    "      {name:<17} {} ({}), worst {}, {} of {} over 1 mm, {:.0} us",
                    fmt(r.rms),
                    fmt(r.crb_rms),
                    fmt(r.worst),
                    r.over_1mm,
                    r.frames - 5,
                    r.us_per_frame
                );
            }
        }
    }
}
