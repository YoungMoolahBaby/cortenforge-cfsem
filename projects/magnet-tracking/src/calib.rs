//! `calib`: calibrating a real board in simulation.
//!
//! Each sensor reads g (B + E) + o: a 3 x 3 gain matrix g (gains and small axis misalignment), an
//! offset o, the background E, at a position off its nominal one by d. A calibration magnet is set
//! at known poses in a jig over the board; for each sensor the 15 unknowns (g, c = g E + o, d) are
//! fitted alone. Two degeneracies make that possible: a gain common to every sensor is the same
//! as a different magnet moment, which the tracker fits anyway, so the calibration magnet's moment
//! can be fixed; and the background at calibration folds into c, while the tracker's own fitted
//! background takes any change after. Readings come from the exact cylinder (Derby and Olbert),
//! so the magnet's shape error is in them; the jig can be off in position and angle.

use nalgebra::{DMatrix, DVector, Matrix3};
use rayon::prelude::*;

use crate::layout::{SENSORS, Sensor, weights};
use crate::magnet::V3;
use crate::solve::{EARTH, MAG, Rng};
use crate::tracker::{Array, Pose, crb_w, dipole, model, solve_w};

/// One sensor's calibration: reading = g * B_at(nominal + d) + c.
#[derive(Clone, Copy, Debug)]
pub struct Cal {
    pub g: Matrix3<f64>,
    pub c: V3,
    pub d: V3,
}

impl Cal {
    fn ideal() -> Self {
        Cal {
            g: Matrix3::identity(),
            c: V3::zeros(),
            d: V3::zeros(),
        }
    }
    fn params(&self) -> DVector<f64> {
        DVector::from_iterator(
            15,
            self.g
                .iter()
                .copied()
                .chain(self.c.iter().copied())
                .chain(self.d.iter().copied()),
        )
    }
    fn from(p: &DVector<f64>) -> Self {
        Cal {
            g: Matrix3::from_iterator(p.iter().take(9).copied()),
            c: V3::new(p[9], p[10], p[11]),
            d: V3::new(p[12], p[13], p[14]),
        }
    }
}

/// A random board's true errors: gains (sd `sg`), misalignment (sd `sm` rad), offsets (sd `so` T)
/// and positions (sd `sd` m), and the background folded into each offset.
fn true_board(rng: &mut Rng, n: usize, sg: f64, sm: f64, so: f64, sd: f64) -> Vec<Cal> {
    (0..n)
        .map(|_| {
            let mut g = Matrix3::identity();
            for i in 0..3 {
                for j in 0..3 {
                    g[(i, j)] += if i == j { sg } else { sm } * rng.normal();
                }
            }
            let o = V3::new(rng.normal(), rng.normal(), rng.normal()) * so;
            Cal {
                g,
                c: g * EARTH + o,
                d: V3::new(rng.normal(), rng.normal(), rng.normal()) * sd,
            }
        })
        .collect()
}

/// A calibration pose: where the jig puts the magnet (nominal) and its axis.
type JigPose = (V3, V3);

/// The jig: magnet centres on a 7 x 7 grid at 15 mm over the board, at 25 and 45 mm, axis up and
/// along x: 196 poses.
fn jig() -> Vec<JigPose> {
    let mut v = Vec::new();
    for ix in -3..=3 {
        for iy in -3..=3 {
            for z in [0.025, 0.045] {
                for u in [V3::z(), V3::x()] {
                    v.push((V3::new(0.015 * ix as f64, 0.015 * iy as f64, z), u));
                }
            }
        }
    }
    v
}

/// Generic Levenberg-Marquardt, central-difference Jacobian with per-parameter steps.
fn lm(f: &dyn Fn(&DVector<f64>) -> DVector<f64>, x0: DVector<f64>, h: &[f64]) -> DVector<f64> {
    let mut x = x0;
    let mut r = f(&x);
    let mut cost = r.norm_squared();
    let mut lambda = 1e-3;
    for _ in 0..200 {
        let mut j = DMatrix::zeros(r.len(), x.len());
        for i in 0..x.len() {
            let mut xp = x.clone();
            xp[i] += h[i];
            let mut xm = x.clone();
            xm[i] -= h[i];
            j.set_column(i, &((f(&xp) - f(&xm)) / (2.0 * h[i])));
        }
        let jtj = j.transpose() * &j;
        let g = j.transpose() * &r;
        let mut improved = false;
        for _ in 0..30 {
            let mut a = jtj.clone();
            for k in 0..a.nrows() {
                a[(k, k)] += lambda * jtj[(k, k)].max(1e-300);
            }
            let Some(ch) = a.cholesky() else {
                lambda *= 10.0;
                continue;
            };
            let xt = &x - ch.solve(&g);
            let rt = f(&xt);
            let ct = rt.norm_squared();
            if ct < cost {
                let rel = (cost - ct) / cost.max(1e-300);
                x = xt;
                r = rt;
                cost = ct;
                lambda = (lambda / 3.0).max(1e-15);
                improved = true;
                if rel < 1e-15 {
                    return x;
                }
                break;
            }
            lambda *= 4.0;
        }
        if !improved {
            break;
        }
    }
    x
}

/// Fit one sensor (nominal position `s`) from its readings at the jig's nominal poses, skipping
/// saturated readings. The magnet is modelled as a dipole of the nominal moment.
fn calibrate(s: &V3, poses: &[JigPose], meas: &[V3], sensor: &Sensor) -> Cal {
    let m0 = MAG.moment();
    let keep: Vec<usize> = (0..poses.len())
        .filter(|&k| meas[k].iter().all(|v| v.abs() <= sensor.range))
        .collect();
    let resid = |p: &DVector<f64>| {
        let cal = Cal::from(p);
        let mut r = DVector::zeros(3 * keep.len());
        for (i, &k) in keep.iter().enumerate() {
            let (pp, u) = poses[k];
            let pred = cal.g * dipole(&pp, &(u * m0), &(s + cal.d)) + cal.c;
            for a in 0..3 {
                r[3 * i + a] = (pred[a] - meas[k][a])
                    / if a == 2 {
                        sensor.noise_z
                    } else {
                        sensor.noise_xy
                    };
            }
        }
        r
    };
    let h: Vec<f64> = (0..15)
        .map(|i| {
            if i < 9 {
                1e-7
            } else if i < 12 {
                1e-12
            } else {
                1e-8
            }
        })
        .collect();
    Cal::from(&lm(&resid, Cal::ideal().params(), &h))
}

/// Calibrate the board from every sensor's readings at the jig poses, each sensor alone, taking
/// the jig's poses as exact.
fn per_sensor(array: &Array, sensor: &Sensor, nominal: &[JigPose], meas: &[Vec<V3>]) -> Vec<Cal> {
    array
        .sensors
        .iter()
        .enumerate()
        .map(|(k, s)| calibrate(s, nominal, &meas[k], sensor))
        .collect()
}

fn tangents(u: &V3) -> (V3, V3) {
    let a = if u.x.abs() < 0.9 { V3::x() } else { V3::y() };
    let e1 = u.cross(&a).normalize();
    (e1, u.cross(&e1))
}

/// Bundle adjustment: every jig pose and every sensor at once (5 per pose, 15 per sensor).
/// Readings weighted by 1 / noise; the jig's nominal poses enter as measurements with its
/// tolerance (`sj` m in position, `sa` rad in angle), the sensors' nominal positions with 1 mm.
/// That fixes what the readings alone leave free: a move, turn or scaling of board and poses
/// together (a scaling's c^-3 in the fields is a common gain), with the jig's errors averaged
/// over its poses. The moment stays fixed (it trades with a common gain). Levenberg-Marquardt on
/// the dense normal equations, built from each reading's 3 x 20 block.
fn bundle(
    array: &Array,
    sensor: &Sensor,
    nominal: &[JigPose],
    meas: &[Vec<V3>],
    sj: f64,
    sa: f64,
) -> Vec<Cal> {
    let m0 = MAG.moment();
    let (np, ns) = (nominal.len(), array.sensors.len());
    let dim = 5 * np + 15 * ns;
    let sd = 1e-3;
    let sig = |a: usize| {
        if a == 2 {
            sensor.noise_z
        } else {
            sensor.noise_xy
        }
    };
    let ok: Vec<Vec<bool>> = meas
        .iter()
        .map(|mk| {
            mk.iter()
                .map(|v| v.iter().all(|x| x.abs() <= sensor.range))
                .collect()
        })
        .collect();
    let step_pose = |(p, u): &JigPose, d: &[f64]| -> JigPose {
        let (e1, e2) = tangents(u);
        (
            p + V3::new(d[0], d[1], d[2]),
            (u + e1 * d[3] + e2 * d[4]).normalize(),
        )
    };
    let reading = |pose: &JigPose, cal: &Cal, k: usize| {
        cal.g * dipole(&pose.0, &(pose.1 * m0), &(array.sensors[k] + cal.d)) + cal.c
    };
    let cost = |poses: &[JigPose], cals: &[Cal]| -> f64 {
        let mut c = 0.0;
        for k in 0..ns {
            for j in 0..np {
                if ok[k][j] {
                    let r = reading(&poses[j], &cals[k], k) - meas[k][j];
                    c += (0..3).map(|a| (r[a] / sig(a)).powi(2)).sum::<f64>();
                }
            }
            c += (cals[k].d / sd).norm_squared();
        }
        for (pose, nom) in poses.iter().zip(nominal) {
            let (n1, n2) = tangents(&nom.1);
            let du = pose.1 - nom.1;
            c += ((pose.0 - nom.0) / sj).norm_squared()
                + (du.dot(&n1) / sa).powi(2)
                + (du.dot(&n2) / sa).powi(2);
        }
        c
    };
    let mut poses = nominal.to_vec();
    let mut cals = per_sensor(array, sensor, nominal, meas);
    let mut c0 = cost(&poses, &cals);
    let mut lambda = 1e-3;
    let hp = [1e-7, 1e-7, 1e-7, 1e-7, 1e-7];
    let hs: Vec<f64> = (0..15)
        .map(|i| {
            if i < 9 {
                1e-7
            } else if i < 12 {
                1e-12
            } else {
                1e-8
            }
        })
        .collect();
    for _ in 0..40 {
        let mut h = DMatrix::<f64>::zeros(dim, dim);
        let mut g = DVector::<f64>::zeros(dim);
        for k in 0..ns {
            let sb = 5 * np + 15 * k;
            for j in 0..np {
                if !ok[k][j] {
                    continue;
                }
                let r0 = reading(&poses[j], &cals[k], k) - meas[k][j];
                let mut jac = DMatrix::<f64>::zeros(3, 20);
                for i in 0..5 {
                    let mut d = [0.0; 5];
                    d[i] = hp[i];
                    let plus = reading(&step_pose(&poses[j], &d), &cals[k], k);
                    d[i] = -hp[i];
                    let minus = reading(&step_pose(&poses[j], &d), &cals[k], k);
                    let col = (plus - minus) / (2.0 * hp[i]);
                    for a in 0..3 {
                        jac[(a, i)] = col[a] / sig(a);
                    }
                }
                let base = cals[k].params();
                for i in 0..15 {
                    let mut x = base.clone();
                    x[i] += hs[i];
                    let plus = reading(&poses[j], &Cal::from(&x), k);
                    x[i] -= 2.0 * hs[i];
                    let minus = reading(&poses[j], &Cal::from(&x), k);
                    let col = (plus - minus) / (2.0 * hs[i]);
                    for a in 0..3 {
                        jac[(a, 5 + i)] = col[a] / sig(a);
                    }
                }
                let r = DVector::from_iterator(3, (0..3).map(|a| r0[a] / sig(a)));
                let idx: Vec<usize> = (0..5)
                    .map(|i| 5 * j + i)
                    .chain((0..15).map(|i| sb + i))
                    .collect();
                let jtj = jac.transpose() * &jac;
                let jtr = jac.transpose() * r;
                for (a, &ia) in idx.iter().enumerate() {
                    g[ia] += jtr[a];
                    for (b, &ib) in idx.iter().enumerate() {
                        h[(ia, ib)] += jtj[(a, b)];
                    }
                }
            }
            // Prior on the sensor's position.
            for i in 0..3 {
                h[(sb + 12 + i, sb + 12 + i)] += 1.0 / (sd * sd);
                g[sb + 12 + i] += cals[k].d[i] / (sd * sd);
            }
        }
        // Priors from the jig.
        for (j, (pose, nom)) in poses.iter().zip(nominal).enumerate() {
            for i in 0..3 {
                h[(5 * j + i, 5 * j + i)] += 1.0 / (sj * sj);
                g[5 * j + i] += (pose.0[i] - nom.0[i]) / (sj * sj);
            }
            let (e1, e2) = tangents(&pose.1);
            let (n1, n2) = tangents(&nom.1);
            let du = pose.1 - nom.1;
            let jr =
                nalgebra::Matrix2::new(e1.dot(&n1), e2.dot(&n1), e1.dot(&n2), e2.dot(&n2)) / sa;
            let rr = nalgebra::Vector2::new(du.dot(&n1), du.dot(&n2)) / sa;
            let jtj = jr.transpose() * jr;
            let jtr = jr.transpose() * rr;
            for a in 0..2 {
                g[5 * j + 3 + a] += jtr[a];
                for b in 0..2 {
                    h[(5 * j + 3 + a, 5 * j + 3 + b)] += jtj[(a, b)];
                }
            }
        }
        let mut improved = false;
        for _ in 0..20 {
            let mut a = h.clone();
            for i in 0..dim {
                a[(i, i)] += lambda * h[(i, i)].max(1e-300);
            }
            let Some(ch) = a.cholesky() else {
                lambda *= 10.0;
                continue;
            };
            let d = -ch.solve(&g);
            let tp: Vec<JigPose> = poses
                .iter()
                .enumerate()
                .map(|(j, p)| step_pose(p, &d.as_slice()[5 * j..5 * j + 5]))
                .collect();
            let tc: Vec<Cal> = (0..ns)
                .map(|k| {
                    Cal::from(
                        &(cals[k].params()
                            + DVector::from_column_slice(
                                &d.as_slice()[5 * np + 15 * k..5 * np + 15 * k + 15],
                            )),
                    )
                })
                .collect();
            let c1 = cost(&tp, &tc);
            if c1 < c0 {
                let rel = (c0 - c1) / c0;
                poses = tp;
                cals = tc;
                c0 = c1;
                lambda = (lambda / 3.0).max(1e-12);
                improved = true;
                if rel < 1e-12 {
                    return cals;
                }
                break;
            }
            lambda *= 4.0;
        }
        if !improved {
            break;
        }
    }
    cals
}

/// Readings of one sensor (true calibration `cal`, nominal position `s`) for a magnet at (p, u),
/// from the exact cylinder, with noise `noise` per axis (already averaged).
fn read(cal: &Cal, s: &V3, p: &V3, u: &V3, noise: (f64, f64), rng: &mut Rng) -> V3 {
    let b = cal.g * MAG.exact(p, u, &(s + cal.d)) + cal.c;
    V3::new(
        b.x + noise.0 * rng.normal(),
        b.y + noise.0 * rng.normal(),
        b.z + noise.1 * rng.normal(),
    )
}

/// Track test poses with a board (true errors `truth`), seen through calibration `est` (None: the
/// nominal board); no noise. Returns the position errors.
fn track_bias(
    truth: &[Cal],
    est: Option<&[Cal]>,
    array: &Array,
    sensor: &Sensor,
    tests: &[Pose],
) -> Vec<f64> {
    let mut rng = Rng(1);
    let seen = Array {
        sensors: array
            .sensors
            .iter()
            .enumerate()
            .map(|(k, s)| s + est.map(|e| e[k].d).unwrap_or_else(V3::zeros))
            .collect(),
    };
    tests
        .iter()
        .map(|t| {
            let mut meas = DVector::zeros(3 * array.sensors.len());
            for (k, s) in array.sensors.iter().enumerate() {
                let raw = read(&truth[k], s, &t.p, &t.u, (0.0, 0.0), &mut rng);
                let b = match est {
                    Some(e) => e[k].g.try_inverse().expect("gain") * (raw - e[k].c),
                    None => raw,
                };
                for a in 0..3 {
                    meas[3 * k + a] = b[a];
                }
            }
            let (w, _) = weights(sensor, &meas);
            let init = Pose {
                g: if est.is_some() { V3::zeros() } else { EARTH },
                ..*t
            };
            (solve_w(&meas, &seen, &init, true, &w).pose.p - t.p).norm()
        })
        .collect()
}

pub fn run() {
    let sensor = SENSORS[0];
    let array = Array::grid(4, 4, 0.03);
    let poses = jig();
    let m0 = MAG.moment();
    println!(
        "Calibration: {}, 4 x 4 at 30 mm; the jig sets the magnet (6 x 5.2 mm) at {} poses (7 x 7 at 15 mm, 25 and 45 mm \
         up, axis up or along x); 15 unknowns a sensor",
        sensor.name,
        poses.len()
    );

    // 1. Exact recovery: no noise, a perfect jig, readings from the dipole itself.
    let mut rng = Rng(0x0000_ddb1_a5e5_u64);
    let board = true_board(&mut rng, 16, 0.01, 0.0175, 0.5e-6, 0.5e-3);
    let mut worst = (0.0_f64, 0.0_f64, 0.0_f64);
    for (k, s) in array.sensors.iter().enumerate() {
        let t = &board[k];
        let meas: Vec<V3> = poses
            .iter()
            .map(|(p, u)| t.g * dipole(p, &(u * m0), &(s + t.d)) + t.c)
            .collect();
        let e = calibrate(s, &poses, &meas, &sensor);
        worst.0 = worst.0.max((e.g - t.g).abs().max());
        worst.1 = worst.1.max((e.c - t.c).norm());
        worst.2 = worst.2.max((e.d - t.d).norm());
    }
    println!(
        "1. No noise, perfect jig, dipole readings, a board with 1% gains, 1 deg misalignment, 0.5 uT offsets, 0.5 mm positions: \
         worst error over 16 sensors: gain {:.1e}, offset {:.1e} T, position {:.1e} m",
        worst.0, worst.1, worst.2
    );

    // 2. Real calibration: noisy readings (16 averaged), the exact cylinder, a jig off by s_jig.
    let tests: Vec<Pose> = [0.02, 0.06, 0.12]
        .iter()
        .flat_map(|&z| {
            [(0.005, -0.003), (0.03, 0.02), (-0.025, 0.035)].map(|(x, y)| Pose {
                p: V3::new(x, y, z),
                u: V3::new(0.3, 0.1, 1.0).normalize(),
                m: m0,
                g: EARTH,
            })
        })
        .collect();
    let floor: Vec<f64> = tests
        .iter()
        .map(|t| {
            let (w, _) = weights(&sensor, &model(t, &array));
            let c = crb_w(t, &array, true, &w).expect("usable");
            (c[(0, 0)] + c[(1, 1)] + c[(2, 2)]).sqrt()
        })
        .collect();
    let noise = (sensor.noise_xy / 4.0, sensor.noise_z / 4.0);
    println!(
        "2. 40 random boards (gains 1%, misalignment 1 deg, offsets 0.5 uT, positions 0.5 mm; assumed for a hand-built board), \
         calibration readings averaged 16x; tracking 9 test poses noise-free, position bias median / 95th percentile, \
         by height (the noise floor's bound in brackets):"
    );
    let heights = [0.02, 0.06, 0.12];
    let floor_at = |z: f64| {
        let v: Vec<f64> = tests
            .iter()
            .zip(&floor)
            .filter(|(t, _)| (t.p.z - z).abs() < 1e-9)
            .map(|(_, f)| *f)
            .collect();
        v.iter().sum::<f64>() / v.len() as f64
    };
    for (label, sj, sa, refit) in [
        ("uncalibrated", f64::NAN, 0.0_f64, false),
        ("perfect jig", 0.0, 0.0, false),
        ("jig 0.05 mm only", 0.05e-3, 0.0, false),
        ("jig 0.2 deg only", 0.0, 0.2, false),
        ("jig 0.05 mm, 0.2 deg", 0.05e-3, 0.2, false),
        ("jig 0.2 mm, 0.5 deg", 0.2e-3, 0.5, false),
        ("same jig, joint solve", 0.2e-3, 0.5, true),
        ("jig 0.5 mm, 1 deg, joint", 0.5e-3, 1.0, true),
    ] {
        let results: Vec<(f64, f64)> = (0..40u64)
            .into_par_iter()
            .flat_map(|b| {
                let mut rng = Rng(0xb0a7_d5ee_d000 ^ (b + 1).wrapping_mul(0x2545_f491_4f6c_dd1d));
                let truth = true_board(&mut rng, 16, 0.01, 0.0175, 0.5e-6, 0.5e-3);
                let est: Option<Vec<Cal>> = if sj.is_nan() {
                    None
                } else {
                    // Where the jig really puts the magnet, against where it says.
                    let real: Vec<JigPose> = poses
                        .iter()
                        .map(|(p, u)| {
                            let tilt =
                                V3::new(rng.normal(), rng.normal(), rng.normal()) * sa.to_radians();
                            (
                                p + V3::new(rng.normal(), rng.normal(), rng.normal()) * sj,
                                (u + tilt).normalize(),
                            )
                        })
                        .collect();
                    let meas: Vec<Vec<V3>> = array
                        .sensors
                        .iter()
                        .enumerate()
                        .map(|(k, s)| {
                            real.iter()
                                .map(|(p, u)| read(&truth[k], s, p, u, noise, &mut rng))
                                .collect()
                        })
                        .collect();
                    Some(if refit {
                        bundle(&array, &sensor, &poses, &meas, sj, sa.to_radians())
                    } else {
                        per_sensor(&array, &sensor, &poses, &meas)
                    })
                };
                let errs = track_bias(&truth, est.as_deref(), &array, &sensor, &tests);
                tests
                    .iter()
                    .zip(errs)
                    .map(|(t, e)| (t.p.z, e))
                    .collect::<Vec<_>>()
            })
            .collect();
        let row: Vec<String> = heights
            .iter()
            .map(|&z| {
                let mut v: Vec<f64> = results
                    .iter()
                    .filter(|r| (r.0 - z).abs() < 1e-9)
                    .map(|r| r.1)
                    .collect();
                v.sort_by(|a, b| a.total_cmp(b));
                let fmt = |x: f64| {
                    if x < 1e-3 {
                        format!("{:.1} um", x * 1e6)
                    } else {
                        format!("{:.2} mm", x * 1e3)
                    }
                };
                format!(
                    "{:.0} mm: {} / {} [{}]",
                    z * 1e3,
                    fmt(v[v.len() / 2]),
                    fmt(v[(0.95 * (v.len() - 1) as f64) as usize]),
                    fmt(floor_at(z))
                )
            })
            .collect();
        println!("   {label:<27} {}", row.join("; "));
    }
}
