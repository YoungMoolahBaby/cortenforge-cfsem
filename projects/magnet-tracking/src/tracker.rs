//! The pose solver: a magnet's position, axis and moment, and a uniform background field, from a
//! magnetometer array's readings, by Levenberg-Marquardt on cfsem's point dipole.

use nalgebra::{DMatrix, DVector};

use crate::magnet::V3;

/// Sensor positions (each reads all three axes, axes aligned with the world).
#[derive(Clone)]
pub struct Array {
    pub sensors: Vec<V3>,
}

impl Array {
    /// nx x ny sensors in the plane z = 0, centred on the origin.
    pub fn grid(nx: usize, ny: usize, pitch: f64) -> Self {
        let mut s = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                s.push(V3::new(
                    (i as f64 - 0.5 * (nx - 1) as f64) * pitch,
                    (j as f64 - 0.5 * (ny - 1) as f64) * pitch,
                    0.0,
                ));
            }
        }
        Array { sensors: s }
    }
}

/// Magnet position, unit axis, moment (A m^2) and the background field (T).
#[derive(Clone, Copy, Debug)]
pub struct Pose {
    pub p: V3,
    pub u: V3,
    pub m: f64,
    pub g: V3,
}

/// The readings a pose gives on an array, as a point dipole plus the background (3 per sensor).
pub fn model(pose: &Pose, array: &Array) -> DVector<f64> {
    let mut out = DVector::zeros(3 * array.sensors.len());
    for (k, s) in array.sensors.iter().enumerate() {
        let b = dipole(&pose.p, &(pose.u * pose.m), s) + pose.g;
        out[3 * k] = b.x;
        out[3 * k + 1] = b.y;
        out[3 * k + 2] = b.z;
    }
    out
}

/// cfsem's point dipole for a moment vector.
pub fn dipole(p: &V3, m: &V3, obs: &V3) -> V3 {
    let b = cfsem::physics::point_source::dipole::flux_density_dipole_scalar(
        (p.x, p.y, p.z),
        (m.x, m.y, m.z),
        0.0,
        (obs.x, obs.y, obs.z),
    );
    V3::new(b.0, b.1, b.2)
}

/// Two unit vectors perpendicular to `u` (and to each other).
fn tangents(u: &V3) -> (V3, V3) {
    let a = if u.x.abs() < 0.9 { V3::x() } else { V3::y() };
    let e1 = u.cross(&a).normalize();
    (e1, u.cross(&e1))
}

/// The pose moved by a parameter step: position (3), axis by two tangent angles (rad), then the
/// moment if `fit_m`, then the background (3).
fn step(pose: &Pose, d: &DVector<f64>, fit_m: bool) -> Pose {
    let (e1, e2) = tangents(&pose.u);
    let o = if fit_m { 6 } else { 5 };
    Pose {
        p: pose.p + V3::new(d[0], d[1], d[2]),
        u: (pose.u + e1 * d[3] + e2 * d[4]).normalize(),
        m: if fit_m { pose.m + d[5] } else { pose.m },
        g: pose.g + V3::new(d[o], d[o + 1], d[o + 2]),
    }
}

pub fn n_params(fit_m: bool) -> usize {
    if fit_m { 9 } else { 8 }
}

/// Jacobian of the readings in those parameters, by central differences.
pub fn jacobian(pose: &Pose, array: &Array, fit_m: bool) -> DMatrix<f64> {
    let np = n_params(fit_m);
    let h: Vec<f64> = (0..np)
        .map(|i| match i {
            0..=2 => 1e-7,
            3 | 4 => 1e-7,
            5 if fit_m => 1e-7 * pose.m.abs().max(1e-6),
            _ => 1e-10,
        })
        .collect();
    let mut j = DMatrix::zeros(3 * array.sensors.len(), np);
    for i in 0..np {
        let mut d = DVector::zeros(np);
        d[i] = h[i];
        let plus = model(&step(pose, &d, fit_m), array);
        d[i] = -h[i];
        let minus = model(&step(pose, &d, fit_m), array);
        j.set_column(i, &((plus - minus) / (2.0 * h[i])));
    }
    j
}

pub struct Fit {
    pub pose: Pose,
    pub cost: f64,
    pub iters: usize,
}

/// Levenberg-Marquardt (Marquardt's diagonal scaling) from `init`. A dipole's field is the same
/// for (u, m) and (-u, -m), so the result is turned to a positive moment.
pub fn solve(meas: &DVector<f64>, array: &Array, init: &Pose, fit_m: bool) -> Fit {
    let mut f = lm(meas, array, init, fit_m);
    if f.pose.m < 0.0 {
        f.pose.m = -f.pose.m;
        f.pose.u = -f.pose.u;
    }
    f
}

fn lm(meas: &DVector<f64>, array: &Array, init: &Pose, fit_m: bool) -> Fit {
    let mut pose = *init;
    let mut r = model(&pose, array) - meas;
    let mut cost = r.norm_squared();
    let mut lambda = 1e-3;
    let mut iters = 0;
    for it in 0..200 {
        iters = it + 1;
        let j = jacobian(&pose, array, fit_m);
        let jtj = j.transpose() * &j;
        let g = j.transpose() * &r;
        let mut improved = false;
        for _ in 0..30 {
            let mut a = jtj.clone();
            for k in 0..a.nrows() {
                a[(k, k)] += lambda * jtj[(k, k)].max(1e-30);
            }
            let Some(ch) = a.cholesky() else {
                lambda *= 10.0;
                continue;
            };
            let d = -ch.solve(&g);
            let trial = step(&pose, &d, fit_m);
            let rt = model(&trial, array) - meas;
            let ct = rt.norm_squared();
            if ct < cost {
                let rel = (cost - ct) / cost.max(1e-300);
                pose = trial;
                r = rt;
                cost = ct;
                lambda = (lambda / 3.0).max(1e-12);
                improved = true;
                if rel < 1e-14 {
                    return Fit { pose, cost, iters };
                }
                break;
            }
            lambda *= 4.0;
        }
        if !improved {
            break;
        }
    }
    Fit { pose, cost, iters }
}

/// Cramer-Rao bound: the parameter covariance sigma^2 (J^T J)^-1 at the true pose.
pub fn crb(pose: &Pose, array: &Array, sigma: f64, fit_m: bool) -> DMatrix<f64> {
    let j = jacobian(pose, array, fit_m);
    (j.transpose() * j).try_inverse().expect("identifiable") * (sigma * sigma)
}
