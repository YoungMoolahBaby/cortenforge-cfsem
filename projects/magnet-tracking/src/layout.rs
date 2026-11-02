//! `layout`: sensor layouts with real magnetometers. Each sensor's noise per axis and its range
//! come from its datasheet; a sensor whose reading on any axis passes its range is dropped from
//! the fit (saturated). Accuracy is the Cramer-Rao bound (moment and background fitted), which
//! `solve` showed the solver reaches; that is checked again here with the weights and drops.

use nalgebra::DVector;
use rayon::prelude::*;

use crate::magnet::V3;
use crate::solve::{EARTH, MAG, Rng};
use crate::tracker::{Array, Pose, crb_w, model, solve_w};

/// A magnetometer: noise per axis (T rms; x and y, then z) and full-scale range (T).
#[derive(Clone, Copy)]
pub struct Sensor {
    pub name: &'static str,
    pub noise_xy: f64,
    pub noise_z: f64,
    pub range: f64,
    pub note: &'static str,
}

/// MEMSIC MMC5983MA (Rev. A, 2019): +-8 G, 0.4 mG "total RMS noise" at BW = 00 (50 Hz output);
/// taken per axis (the pessimistic reading). ST LIS3MDL (DocID024204 Rev 6): RMS noise 3.2 mG on x
/// and y, 4.1 mG on z, ultra-high-performance mode at +-12 G.
pub const SENSORS: [Sensor; 2] = [
    Sensor {
        name: "MMC5983MA",
        noise_xy: 0.04e-6,
        noise_z: 0.04e-6,
        range: 800e-6,
        note: "+-8 G, 50 Hz",
    },
    Sensor {
        name: "LIS3MDL",
        noise_xy: 0.32e-6,
        noise_z: 0.41e-6,
        range: 1200e-6,
        note: "+-12 G, UHP",
    },
];

/// Weights (1 / noise) for a pose's readings, 0 for every axis of a saturated sensor; and how many
/// sensors are left.
pub fn weights(sensor: &Sensor, clean: &DVector<f64>) -> (DVector<f64>, usize) {
    let n = clean.len() / 3;
    let mut w = DVector::zeros(clean.len());
    let mut kept = 0;
    for k in 0..n {
        if (0..3).all(|a| clean[3 * k + a].abs() <= sensor.range) {
            kept += 1;
            w[3 * k] = 1.0 / sensor.noise_xy;
            w[3 * k + 1] = 1.0 / sensor.noise_xy;
            w[3 * k + 2] = 1.0 / sensor.noise_z;
        }
    }
    (w, kept)
}

/// Position error bound (m; root of the trace of the position covariance), or None if unusable.
fn pos_bound(sensor: &Sensor, array: &Array, pose: &Pose) -> Option<f64> {
    let clean = model(pose, array);
    let (w, kept) = weights(sensor, &clean);
    if kept < 4 {
        return None;
    }
    let c = crb_w(pose, array, true, &w)?;
    Some((c[(0, 0)] + c[(1, 1)] + c[(2, 2)]).sqrt())
}

pub fn run() {
    let m0 = MAG.moment();
    let heights = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12];
    let axes = [V3::z(), V3::x(), V3::new(1.0, 1.0, 1.0).normalize()];
    println!(
        "Sensors: {}",
        SENSORS
            .iter()
            .map(|s| format!(
                "{} ({}): {:.2} / {:.2} uT rms (xy / z), range {:.0} uT",
                s.name,
                s.note,
                s.noise_xy * 1e6,
                s.noise_z * 1e6,
                s.range * 1e6
            ))
            .collect::<Vec<_>>()
            .join("; ")
    );
    println!(
        "Magnet as in `solve` (6 x 5.2 mm, {m0:.4} A m^2); moment and background fitted; saturated sensors dropped"
    );

    // 1. The weighted solver against the weighted bound, with drops.
    println!("1. Weighted solver vs the weighted Cramer-Rao bound, 4000 trials, 4 x 4 at 30 mm:");
    let array = Array::grid(4, 4, 0.03);
    for (sensor, p, u) in [
        (SENSORS[0], V3::new(0.015, 0.015, 0.018), V3::z()),
        (
            SENSORS[1],
            V3::new(0.01, -0.01, 0.05),
            V3::new(1.0, 0.3, 0.6).normalize(),
        ),
    ] {
        let truth = Pose {
            p,
            u,
            m: m0,
            g: EARTH,
        };
        let clean = model(&truth, &array);
        let (w, kept) = weights(&sensor, &clean);
        let c = crb_w(&truth, &array, true, &w).expect("usable");
        let est: Vec<V3> = (0..4000u64)
            .into_par_iter()
            .map(|k| {
                let mut r =
                    Rng(0x9e37_79b9_7f4a_7c15 ^ (k + 7).wrapping_mul(0x2545_f491_4f6c_dd1d));
                let meas = DVector::from_iterator(
                    clean.len(),
                    clean.iter().enumerate().map(|(i, v)| {
                        v + if i % 3 == 2 {
                            sensor.noise_z
                        } else {
                            sensor.noise_xy
                        } * r.normal()
                    }),
                );
                solve_w(&meas, &array, &truth, true, &w).pose.p - truth.p
            })
            .collect();
        let n = est.len() as f64;
        let mean = est.iter().fold(V3::zeros(), |a, e| a + e) / n;
        let sd = |i: usize| {
            (est.iter().map(|e| (e[i] - mean[i]).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()
        };
        println!(
            "   {:<9} magnet {:.0} mm up, {kept} of 16 sensors unsaturated: x {:.2} / {:.2} um, y {:.2} / {:.2} um, z {:.2} / {:.2} um (spread / bound)",
            sensor.name,
            p.z * 1e3,
            sd(0) * 1e6,
            c[(0, 0)].sqrt() * 1e6,
            sd(1) * 1e6,
            c[(1, 1)].sqrt() * 1e6,
            sd(2) * 1e6,
            c[(2, 2)].sqrt() * 1e6
        );
    }

    // 2. Layouts over the workspace.
    println!(
        "2. Position error bound over the workspace (magnet x, y in -40..40 mm in 20 mm steps, axis z, x or diagonal), \
         95th percentile per height (unusable poses count as failures); reach = heights with it under 1 mm"
    );
    for sensor in SENSORS {
        println!("   {} ({}):", sensor.name, sensor.note);
        println!(
            "      layout          board   20 mm      40 mm      60 mm      80 mm      100 mm     120 mm     within 1 mm  most saturated"
        );
        for n in [3usize, 4, 5, 6] {
            for pitch in [0.02, 0.03, 0.04, 0.05] {
                let array = Array::grid(n, n, pitch);
                let mut cells = Vec::new();
                let mut ok_heights = Vec::new();
                let mut most_dropped = 0;
                for &z in &heights {
                    let mut errs: Vec<f64> = Vec::new();
                    for ix in -2..=2 {
                        for iy in -2..=2 {
                            for u in &axes {
                                let pose = Pose {
                                    p: V3::new(0.02 * ix as f64, 0.02 * iy as f64, z),
                                    u: *u,
                                    m: m0,
                                    g: EARTH,
                                };
                                errs.push(
                                    pos_bound(&sensor, &array, &pose).unwrap_or(f64::INFINITY),
                                );
                                let (_, kept) = weights(&sensor, &model(&pose, &array));
                                most_dropped = most_dropped.max(n * n - kept);
                            }
                        }
                    }
                    errs.sort_by(|a, b| a.total_cmp(b));
                    let p95 = errs[(0.95 * (errs.len() - 1) as f64).round() as usize];
                    let unusable = errs.iter().filter(|e| e.is_infinite()).count();
                    cells.push(if p95.is_infinite() {
                        format!("{:<10}", format!("sat {}/{}", unusable, errs.len()))
                    } else if p95 < 1e-3 {
                        format!("{:<10}", format!("{:.0} um", p95 * 1e6))
                    } else {
                        format!("{:<10}", format!("{:.1} mm", p95 * 1e3))
                    });
                    if p95 < 1e-3 {
                        ok_heights.push(z);
                    }
                }
                let reach = match (ok_heights.first(), ok_heights.last()) {
                    (Some(a), Some(b)) => format!("{:.0}-{:.0} mm", a * 1e3, b * 1e3),
                    _ => "none".into(),
                };
                println!(
                    "      {n} x {n} at {:.0} mm  {:>3.0} mm  {} {reach:<12} {most_dropped} of {}",
                    pitch * 1e3,
                    (n - 1) as f64 * pitch * 1e3,
                    cells.join(" "),
                    n * n
                );
            }
        }
    }
    calibration();
}

/// Section 3, calibration errors on a real board: each sensor axis with a random gain error and each
/// sensor sitting a random distance from where the solver thinks it is; no noise, so the position
/// error is all bias. MMC5983MA, 4 x 4 at 30 mm, 300 random boards.
pub fn calibration() {
    let m0 = MAG.moment();
    let sensor = SENSORS[0];
    let array = Array::grid(4, 4, 0.03);
    println!(
        "3. Calibration: {}, 4 x 4 at 30 mm, 300 random boards, no noise: position bias (median / 95th percentile) \
         against the noise's bound",
        sensor.name
    );
    for (label, sg, sp) in [
        ("gain 0.1% per axis", 1e-3, 0.0),
        ("gain 1% per axis", 1e-2, 0.0),
        ("sensors 0.1 mm off", 0.0, 1e-4),
        ("sensors 0.5 mm off", 0.0, 5e-4),
    ] {
        let mut row = Vec::new();
        for z in [0.02, 0.06, 0.12] {
            let truth = Pose {
                p: V3::new(0.005, -0.003, z),
                u: V3::new(0.3, 0.0, 1.0).normalize(),
                m: m0,
                g: EARTH,
            };
            let (w, _) = weights(&sensor, &model(&truth, &array));
            let c = crb_w(&truth, &array, true, &w).expect("usable");
            let bound = (c[(0, 0)] + c[(1, 1)] + c[(2, 2)]).sqrt();
            let mut errs: Vec<f64> = (0..300u64)
                .into_par_iter()
                .map(|b| {
                    let mut r =
                        Rng(0x51_7cc1_b727_220a ^ (b + 3).wrapping_mul(0x2545_f491_4f6c_dd1d));
                    let mut meas = DVector::zeros(3 * array.sensors.len());
                    for (k, s) in array.sensors.iter().enumerate() {
                        let at = s + V3::new(r.normal(), r.normal(), r.normal()) * sp;
                        let bfield =
                            crate::tracker::dipole(&truth.p, &(truth.u * truth.m), &at) + truth.g;
                        for a in 0..3 {
                            meas[3 * k + a] = bfield[a] * (1.0 + sg * r.normal());
                        }
                    }
                    (solve_w(&meas, &array, &truth, true, &w).pose.p - truth.p).norm()
                })
                .collect();
            errs.sort_by(|a, b| a.total_cmp(b));
            let (med, p95) = (
                errs[errs.len() / 2],
                errs[(0.95 * (errs.len() - 1) as f64) as usize],
            );
            row.push(format!(
                "{:.0} mm: {:.0} / {:.0} um ({:.1}x)",
                z * 1e3,
                med * 1e6,
                p95 * 1e6,
                med / bound
            ));
        }
        println!("   {label:<20} {}", row.join("; "));
    }
}
