//! The method of regularized Stokeslets (Cortez, Fauci and Medovikov, Phys. Fluids 17, 031504
//! (2005)): the reference model for a rigid body in Stokes flow.
//!
//! A force f spread by the blob phi(r) = 15 eps^4 / (8 pi (r^2 + eps^2)^(7/2)) moves the fluid at
//! u(x) = [f (r^2 + 2 eps^2) + (f.d) d] / (8 pi mu (r^2 + eps^2)^(3/2)), d = x - x0. Points on the
//! body's surface carry unknown forces; requiring the body's rigid velocity at every point gives a
//! symmetric positive definite system, solved here by conjugate gradients, matrix-free, for all
//! right-hand sides at once.

use std::f64::consts::PI;

use nalgebra::Matrix6;
use rayon::prelude::*;

use crate::helix::V3;

/// eps = EPS_C h, h the spacing between points: the c with the least worst sphere error at every
/// spacing (`check` section 2, which also confirms this constant).
pub const EPS_C: f64 = 0.3;

/// The model's error falls about as the spacing h (`check` section 3), so from 12 and 16 points
/// round a wire (h in the ratio 4 : 3) the fine-spacing limit is R = 4 R16 - 3 R12.
pub fn extrapolate(r12: &Matrix6<f64>, r16: &Matrix6<f64>) -> Matrix6<f64> {
    extrapolate_n(r12, 12, r16, 16)
}

/// The same for any two resolutions: R = (n_hi R_hi - n_lo R_lo) / (n_hi - n_lo).
pub fn extrapolate_n(
    r_lo: &Matrix6<f64>,
    n_lo: usize,
    r_hi: &Matrix6<f64>,
    n_hi: usize,
) -> Matrix6<f64> {
    (r_hi * n_hi as f64 - r_lo * n_lo as f64) / (n_hi - n_lo) as f64
}

/// Velocity at `x` from a regularized force `f` at `x0`.
pub fn velocity(x: &V3, x0: &V3, f: &V3, eps: f64, mu: f64) -> V3 {
    let d = x - x0;
    let r2 = d.norm_squared();
    let big = r2 + eps * eps;
    let inv = 1.0 / (8.0 * PI * mu * big * big.sqrt());
    f * ((r2 + 2.0 * eps * eps) * inv) + d * (inv * d.dot(f))
}

/// u_i = sum_j S(x_i - x_j) f_j for several force sets at once (one kernel evaluation a pair).
fn apply(pts: &[V3], eps: f64, mu: f64, fs: &[Vec<V3>]) -> Vec<Vec<V3>> {
    let k = fs.len();
    let e2 = eps * eps;
    let c = 1.0 / (8.0 * PI * mu);
    let rows: Vec<Vec<V3>> = pts
        .par_iter()
        .map(|xi| {
            let mut acc = vec![V3::zeros(); k];
            for (j, xj) in pts.iter().enumerate() {
                let d = xi - xj;
                let r2 = d.norm_squared();
                let big = r2 + e2;
                let inv = c / (big * big.sqrt());
                let a = (r2 + 2.0 * e2) * inv;
                for (q, f) in fs.iter().enumerate() {
                    let fj = &f[j];
                    acc[q] += fj * a + d * (inv * d.dot(fj));
                }
            }
            acc
        })
        .collect();
    (0..k)
        .map(|q| rows.iter().map(|r| r[q]).collect())
        .collect()
}

fn dot(a: &[V3], b: &[V3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x.dot(y)).sum()
}

/// Forces on the points that give the velocities `rhs` (several sets), by conjugate gradients to
/// a relative residual `tol`. Returns the forces and the iterations taken.
pub fn solve(pts: &[V3], eps: f64, mu: f64, rhs: &[Vec<V3>], tol: f64) -> (Vec<Vec<V3>>, usize) {
    let k = rhs.len();
    let n = pts.len();
    let mut x = vec![vec![V3::zeros(); n]; k];
    let mut r: Vec<Vec<V3>> = rhs.to_vec();
    let mut p = r.clone();
    let mut rr: Vec<f64> = r.iter().map(|v| dot(v, v)).collect();
    let bn: Vec<f64> = rr.iter().map(|v| v.sqrt()).collect();
    let mut it = 0;
    while it < 20_000 {
        let done: Vec<bool> = (0..k).map(|q| rr[q].sqrt() <= tol * bn[q]).collect();
        if done.iter().all(|&d| d) {
            break;
        }
        let ap = apply(pts, eps, mu, &p);
        for q in 0..k {
            if done[q] {
                continue;
            }
            let alpha = rr[q] / dot(&p[q], &ap[q]);
            for i in 0..n {
                x[q][i] += p[q][i] * alpha;
                r[q][i] -= ap[q][i] * alpha;
            }
            let rr_new = dot(&r[q], &r[q]);
            let beta = rr_new / rr[q];
            for i in 0..n {
                p[q][i] = r[q][i] + p[q][i] * beta;
            }
            rr[q] = rr_new;
        }
        it += 1;
    }
    (x, it)
}

/// The 6 x 6 resistance matrix about `centre`: [F; T] = R [U; Omega], F and T the force and torque
/// the body exerts on the fluid. Returns R and the iterations taken.
pub fn resistance(pts: &[V3], eps: f64, mu: f64, centre: &V3, tol: f64) -> (Matrix6<f64>, usize) {
    let rhs: Vec<Vec<V3>> = (0..6)
        .map(|col| {
            pts.iter()
                .map(|x| {
                    let mut e = V3::zeros();
                    e[col % 3] = 1.0;
                    if col < 3 { e } else { e.cross(&(x - centre)) }
                })
                .collect()
        })
        .collect();
    let (f, it) = solve(pts, eps, mu, &rhs, tol);
    let mut r = Matrix6::zeros();
    for col in 0..6 {
        let mut force = V3::zeros();
        let mut torque = V3::zeros();
        for (x, fi) in pts.iter().zip(&f[col]) {
            force += fi;
            torque += (x - centre).cross(fi);
        }
        for i in 0..3 {
            r[(i, col)] = force[i];
            r[(3 + i, col)] = torque[i];
        }
    }
    (r, it)
}
