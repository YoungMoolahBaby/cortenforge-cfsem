//! `solve`: the pose solver on a simulated array: exact recovery without noise, then noise against
//! the Cramer-Rao bound, then the bias a real (non-dipole) magnet puts into a dipole solver.

use std::f64::consts::PI;
use std::time::Instant;

use nalgebra::DVector;
use rayon::prelude::*;

use crate::magnet::{Magnet, V3};
use crate::tracker::{Array, Pose, crb, model, solve};

/// The tracking magnet: 6 mm across, 5.2 mm long (L/D 0.867, the dipole-like aspect), N42.
pub const MAG: Magnet = Magnet {
    radius: 3e-3,
    length: 5.2e-3,
    br: 1.3,
};
/// Placeholder noise per axis (T): every error below scales with it.
pub const SIGMA: f64 = 1e-6;
/// A background field like the Earth's (T).
pub const EARTH: V3 = V3::new(20e-6, -5e-6, -40e-6);

/// xorshift64 uniform in [0, 1) and Box-Muller normals.
pub struct Rng(pub u64);
impl Rng {
    pub fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn normal(&mut self) -> f64 {
        let (a, b) = (self.uniform().max(1e-300), self.uniform());
        (-2.0 * a.ln()).sqrt() * (2.0 * PI * b).cos()
    }
    pub fn direction(&mut self) -> V3 {
        V3::new(self.normal(), self.normal(), self.normal()).normalize()
    }
}

fn angle(a: &V3, b: &V3) -> f64 {
    a.normalize().dot(&b.normalize()).clamp(-1.0, 1.0).acos()
}

/// Cold start: the best of a 3 x 3 x 3 grid of positions and six axis directions.
pub fn acquire(meas: &DVector<f64>, array: &Array, m0: f64) -> crate::tracker::Fit {
    let mut starts = Vec::new();
    for x in [-0.03, 0.0, 0.03] {
        for y in [-0.03, 0.0, 0.03] {
            for z in [0.025, 0.05, 0.075] {
                for u in [V3::x(), -V3::x(), V3::y(), -V3::y(), V3::z(), -V3::z()] {
                    starts.push(Pose {
                        p: V3::new(x, y, z),
                        u,
                        m: m0,
                        g: V3::zeros(),
                    });
                }
            }
        }
    }
    starts
        .par_iter()
        .map(|s| solve(meas, array, s, true))
        .reduce_with(|a, b| if b.cost < a.cost { b } else { a })
        .expect("starts")
}

pub fn run() {
    let array = Array::grid(4, 4, 0.03);
    let m0 = MAG.moment();
    println!(
        "Magnet {:.0} x {:.1} mm (L/D {:.3}), Br {} T, moment {:.4} A m^2; array 4 x 4 at 30 mm pitch (90 x 90 mm), \
         noise {:.1} uT per axis (placeholder), background ({:.0}, {:.0}, {:.0}) uT; 9 unknowns (position, axis, moment, background)",
        2e3 * MAG.radius,
        1e3 * MAG.length,
        MAG.length / (2.0 * MAG.radius),
        MAG.br,
        m0,
        SIGMA * 1e6,
        EARTH.x * 1e6,
        EARTH.y * 1e6,
        EARTH.z * 1e6
    );

    // 1. Noise-free recovery, readings from the dipole itself.
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let truths: Vec<Pose> = (0..50)
        .map(|_| Pose {
            p: V3::new(
                0.08 * (rng.uniform() - 0.5),
                0.08 * (rng.uniform() - 0.5),
                0.02 + 0.06 * rng.uniform(),
            ),
            u: rng.direction(),
            m: m0,
            g: EARTH,
        })
        .collect();
    let (mut worst_p, mut worst_a, mut its) = (0.0_f64, 0.0_f64, 0);
    for t in &truths {
        let meas = model(t, &array);
        let off = rng.direction() * 0.01;
        let tilt = (t.u + rng.direction() * (30f64.to_radians())).normalize();
        let init = Pose {
            p: t.p + off,
            u: tilt,
            m: 1.1 * m0,
            g: V3::zeros(),
        };
        let f = solve(&meas, &array, &init, true);
        worst_p = worst_p.max((f.pose.p - t.p).norm());
        worst_a = worst_a.max(angle(&f.pose.u, &t.u));
        its += f.iters;
    }
    println!(
        "1. No noise, 50 poses (x, y within 40 mm of centre, 20-80 mm up, any axis), started 10 mm and ~30 deg off with \
         the moment 10% high and no background: worst position error {:.1e} m, axis {:.1e} rad; {:.1} iterations on average",
        worst_p,
        worst_a,
        its as f64 / 50.0
    );
    let t0 = Instant::now();
    let cold: Vec<(f64, f64)> = truths
        .iter()
        .map(|t| {
            let f = acquire(&model(t, &array), &array, m0);
            ((f.pose.p - t.p).norm(), angle(&f.pose.u, &t.u))
        })
        .collect();
    let found = cold.iter().filter(|c| c.0 < 1e-6 && c.1 < 1e-6).count();
    println!(
        "   cold start (best of 162 starts: 27 positions x 6 axes): {found} of 50 found to 1 um and 1 urad; {:.0} ms an acquisition",
        t0.elapsed().as_secs_f64() * 1e3 / 50.0
    );
    for (t, c) in truths.iter().zip(&cold) {
        if c.0 >= 1e-6 || c.1 >= 1e-6 {
            println!(
                "      missed: magnet at ({:.1}, {:.1}, {:.1}) mm: off by {:.1} mm, {:.1} deg",
                t.p.x * 1e3,
                t.p.y * 1e3,
                t.p.z * 1e3,
                c.0 * 1e3,
                c.1.to_degrees()
            );
        }
    }

    // 2. Noise: Monte Carlo spread of the estimates against the Cramer-Rao bound.
    println!(
        "2. With noise ({:.1} uT per axis), 1000 trials a pose, solved from the truth, vs the Cramer-Rao bound sigma^2 (J^T J)^-1:",
        SIGMA * 1e6
    );
    let cases = [
        (
            "centre, 30 mm up, axis up",
            V3::new(0.0, 0.0, 0.03),
            V3::z(),
        ),
        (
            "centre, 70 mm up, axis up",
            V3::new(0.0, 0.0, 0.07),
            V3::z(),
        ),
        (
            "off the corner, 40 mm up, tilted",
            V3::new(0.04, 0.03, 0.04),
            V3::new(1.0, -0.5, 0.7).normalize(),
        ),
    ];
    for (name, p, u) in cases {
        let truth = Pose {
            p,
            u,
            m: m0,
            g: EARTH,
        };
        let clean = model(&truth, &array);
        let c = crb(&truth, &array, SIGMA, true);
        let c_known = crb(&truth, &array, SIGMA, false);
        let est: Vec<(V3, f64)> = (0..1000u64)
            .into_par_iter()
            .map(|k| {
                let mut r =
                    Rng(0x9e37_79b9_7f4a_7c15 ^ (k + 1).wrapping_mul(0x2545_f491_4f6c_dd1d));
                let meas = DVector::from_iterator(
                    clean.len(),
                    clean.iter().map(|v| v + SIGMA * r.normal()),
                );
                let f = solve(&meas, &array, &truth, true);
                (f.pose.p - truth.p, angle(&f.pose.u, &truth.u))
            })
            .collect();
        let n = est.len() as f64;
        let mean = est.iter().fold(V3::zeros(), |a, e| a + e.0) / n;
        let sd = |i: usize| {
            (est.iter().map(|e| (e.0[i] - mean[i]).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()
        };
        let ang = (est.iter().map(|e| e.1 * e.1).sum::<f64>() / n).sqrt();
        let bound = |i: usize| c[(i, i)].sqrt();
        println!(
            "   {name:<33} x {:.1} / {:.1} um, y {:.1} / {:.1} um, z {:.1} / {:.1} um, axis {:.3} / {:.3} deg (spread / bound); \
             bias {:.2} um; z bound with the moment known {:.1} um",
            sd(0) * 1e6,
            bound(0) * 1e6,
            sd(1) * 1e6,
            bound(1) * 1e6,
            sd(2) * 1e6,
            bound(2) * 1e6,
            ang.to_degrees(),
            (c[(3, 3)] + c[(4, 4)]).sqrt().to_degrees(),
            mean.norm() * 1e6,
            c_known[(2, 2)].sqrt() * 1e6
        );
    }

    // 3. Model error: readings from the real cylinder (closed form), solved as a dipole.
    println!(
        "3. A real magnet read by a dipole solver (no noise): position and axis bias, magnets of the same volume and Br"
    );
    let vol = PI * MAG.radius * MAG.radius * MAG.length;
    for ld in [0.5, 0.867, 1.0, 2.0] {
        let r = (vol / (2.0 * PI * ld)).cbrt();
        let mag = Magnet {
            radius: r,
            length: 2.0 * r * ld,
            br: MAG.br,
        };
        let mut row = Vec::new();
        for z in [0.015, 0.025, 0.04, 0.08] {
            let u = V3::new(0.5, 0.0, 1.0).normalize();
            let p = V3::new(0.01, -0.005, z);
            let meas = DVector::from_iterator(
                3 * array.sensors.len(),
                array.sensors.iter().flat_map(|s| {
                    let b = mag.exact(&p, &u, s) + EARTH;
                    [b.x, b.y, b.z]
                }),
            );
            let truth = Pose {
                p,
                u,
                m: mag.moment(),
                g: EARTH,
            };
            let f = solve(&meas, &array, &truth, true);
            let bound = crb(&truth, &array, SIGMA, true);
            let pos_bound = (bound[(0, 0)] + bound[(1, 1)] + bound[(2, 2)]).sqrt();
            row.push(format!(
                "{:.0} mm: {:.1} um, {:.3} deg ({:.2} x the noise's)",
                z * 1e3,
                (f.pose.p - p).norm() * 1e6,
                angle(&f.pose.u, &u).to_degrees(),
                (f.pose.p - p).norm() / pos_bound
            ));
        }
        println!("   L/D {ld:<5}  {}", row.join("; "));
    }
}
