//! Resistive force theory: each piece of the wire feels a drag set by its own velocity alone,
//! xi_par along the wire and xi_perp across it (Gray and Hancock, J. Exp. Biol. 32, 802 (1955);
//! Lighthill, SIAM Rev. 18, 161 (1976)). Swimming comes from xi_perp > xi_par.
//!
//! Also here: the swimming motion from a 6 x 6 resistance matrix.

use std::f64::consts::PI;

use nalgebra::{Matrix3, Matrix5, Matrix6, Vector3, Vector5, Vector6};

use crate::helix::{Helix, V3};

/// Drag coefficients per unit length (N s / m^2).
#[derive(Clone, Copy, Debug)]
pub struct Coeffs {
    pub perp: f64,
    pub par: f64,
}

/// Gray and Hancock: xi_perp = 4 pi mu / (ln(2 lambda / b) + 1/2), xi_par = 2 pi mu /
/// (ln(2 lambda / b) - 1/2), lambda the wavelength along the axis (the pitch).
pub fn gray_hancock(h: &Helix, mu: f64) -> Coeffs {
    let l = (2.0 * h.pitch / h.wire).ln();
    Coeffs {
        perp: 4.0 * PI * mu / (l + 0.5),
        par: 2.0 * PI * mu / (l - 0.5),
    }
}

/// Lighthill: xi_perp = 4 pi mu / (ln(2 q / b) + 1/2), xi_par = 2 pi mu / ln(2 q / b), q = 0.09
/// Lambda, Lambda the wavelength measured along the wire (its length per turn).
pub fn lighthill(h: &Helix, mu: f64) -> Coeffs {
    let q = 0.09 * 2.0 * PI * h.speed();
    let l = (2.0 * q / h.wire).ln();
    Coeffs {
        perp: 4.0 * PI * mu / (l + 0.5),
        par: 2.0 * PI * mu / l,
    }
}

/// The 6 x 6 resistance matrix about the helix's centre ([F; T] on the fluid = R [U; Omega]),
/// midpoint rule over `n` pieces. Includes the wire spinning about its own centreline (torque
/// 4 pi mu b^2 per length per unit spin), small next to the rest.
pub fn resistance(h: &Helix, k: Coeffs, mu: f64, n: usize) -> Matrix6<f64> {
    let c = h.centre();
    let ds = h.length() / n as f64;
    let mut r = Matrix6::zeros();
    for i in 0..n {
        let t = h.t_max() * (i as f64 + 0.5) / n as f64;
        let x = h.point(t) - c;
        let (tg, _, _) = h.frame(t);
        for col in 0..6 {
            let mut u = V3::zeros();
            let mut w = V3::zeros();
            if col < 3 {
                u[col] = 1.0;
            } else {
                w[col - 3] = 1.0;
            }
            let v = u + w.cross(&x);
            let vt = v.dot(&tg);
            let f = tg * (k.par * vt) + (v - tg * vt) * k.perp;
            let tq = x.cross(&f) + tg * (4.0 * PI * mu * h.wire * h.wire * w.dot(&tg));
            for j in 0..3 {
                r[(j, col)] += f[j] * ds;
                r[(3 + j, col)] += tq[j] * ds;
            }
        }
    }
    r
}

/// Free swimming turned at unit rate about z: no net force, no torque across the axis (the
/// magnet's torque is along z). Returns the velocity, the wobble (Omega_x, Omega_y) and the
/// torque the magnet must supply, all per rad/s.
pub struct Swim {
    pub u: V3,
    pub wobble: (f64, f64),
    pub torque: f64,
}

pub fn swim(r: &Matrix6<f64>) -> Swim {
    let a: Matrix5<f64> = r.fixed_view::<5, 5>(0, 0).into_owned();
    let b: Vector5<f64> = -r.fixed_view::<5, 1>(0, 5).into_owned();
    let y = a.lu().solve(&b).expect("resistance block invertible");
    let full = Vector6::new(y[0], y[1], y[2], y[3], y[4], 1.0);
    Swim {
        u: V3::new(y[0], y[1], y[2]),
        wobble: (y[3], y[4]),
        torque: (r.row(5) * full)[0],
    }
}

/// Advance per turn as a fraction of the pitch (1 would be a screw in wood).
pub fn advance(h: &Helix, s: &Swim) -> f64 {
    s.u.z * 2.0 * PI / h.pitch
}

/// Turned at unit rate about z with the axis held (no wobble): no net force; the torque about
/// x and y is whatever holding the axis takes.
pub fn swim_locked(r: &Matrix6<f64>) -> Swim {
    let a: Matrix3<f64> = r.fixed_view::<3, 3>(0, 0).into_owned();
    let b: Vector3<f64> = -r.fixed_view::<3, 1>(0, 5).into_owned();
    let u = a.lu().solve(&b).expect("translation block invertible");
    Swim {
        u,
        wobble: (0.0, 0.0),
        torque: r[(5, 0)] * u.x + r[(5, 1)] * u.y + r[(5, 2)] * u.z + r[(5, 5)],
    }
}
