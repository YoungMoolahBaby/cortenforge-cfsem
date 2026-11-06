//! Modular stellarator coils as Fourier-series curves, and their field through cfsem's
//! straight-filament Biot-Savart (`flux_density_linear_filament`).
//!
//! A coil centreline is x(t) = sum_m xc_m cos(m t) + xs_m sin(m t) (same for y and z),
//! t in [0, 2 pi), the form stellarator coil optimizers use. Those coefficients are the
//! numbers a later optimizer will change.

use std::f64::consts::PI;

use cfsem::MU_0;
use cfsem::physics::linear_filament::flux_density_linear_filament;

/// A closed curve as a Fourier series in each Cartesian component.
#[derive(Clone, Debug)]
pub struct FourierCurve {
    /// cos coefficients, [component][m], m = 0..=order.
    pub c: [Vec<f64>; 3],
    /// sin coefficients, [component][m], m = 0..=order (m = 0 unused).
    pub s: [Vec<f64>; 3],
}

impl FourierCurve {
    pub fn zeros(order: usize) -> Self {
        let z = vec![0.0; order + 1];
        FourierCurve {
            c: [z.clone(), z.clone(), z.clone()],
            s: [z.clone(), z.clone(), z],
        }
    }

    pub fn order(&self) -> usize {
        self.c[0].len() - 1
    }

    pub fn point(&self, t: f64) -> [f64; 3] {
        let mut p = [0.0; 3];
        for (k, pk) in p.iter_mut().enumerate() {
            *pk = self.c[k][0];
            for m in 1..=self.order() {
                let (s, c) = (m as f64 * t).sin_cos();
                *pk += self.c[k][m] * c + self.s[k][m] * s;
            }
        }
        p
    }

    pub fn tangent(&self, t: f64) -> [f64; 3] {
        let mut d = [0.0; 3];
        for (k, dk) in d.iter_mut().enumerate() {
            for m in 1..=self.order() {
                let mf = m as f64;
                let (s, c) = (mf * t).sin_cos();
                *dk += mf * (-self.c[k][m] * s + self.s[k][m] * c);
            }
        }
        d
    }

    pub fn second(&self, t: f64) -> [f64; 3] {
        let mut d = [0.0; 3];
        for (k, dk) in d.iter_mut().enumerate() {
            for m in 1..=self.order() {
                let mf = m as f64;
                let (s, c) = (mf * t).sin_cos();
                *dk -= mf * mf * (self.c[k][m] * c + self.s[k][m] * s);
            }
        }
        d
    }

    /// Curvature |x' x x''| / |x'|^3 at t (1/m).
    pub fn curvature(&self, t: f64) -> f64 {
        let (d1, d2) = (self.tangent(t), self.second(t));
        let c = [
            d1[1] * d2[2] - d1[2] * d2[1],
            d1[2] * d2[0] - d1[0] * d2[2],
            d1[0] * d2[1] - d1[1] * d2[0],
        ];
        let n1 = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
        (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt() / n1.powi(3)
    }

    /// The same curve padded (or truncated) to `order`.
    pub fn with_order(&self, order: usize) -> Self {
        let mut f = FourierCurve::zeros(order);
        for k in 0..3 {
            for m in 0..=order.min(self.order()) {
                f.c[k][m] = self.c[k][m];
                f.s[k][m] = self.s[k][m];
            }
        }
        f
    }

    /// Rotated by `angle` about the z axis.
    pub fn rotate_z(&self, angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        let mut f = self.clone();
        for m in 0..=self.order() {
            f.c[0][m] = c * self.c[0][m] - s * self.c[1][m];
            f.c[1][m] = s * self.c[0][m] + c * self.c[1][m];
            f.s[0][m] = c * self.s[0][m] - s * self.s[1][m];
            f.s[1][m] = s * self.s[0][m] + c * self.s[1][m];
        }
        f
    }

    /// The stellarator-symmetric partner: rotated by pi about the x axis, (x, y, z) ->
    /// (x, -y, -z), and run backwards (t -> -t) so the partner's current keeps B_phi's sign.
    pub fn mirror(&self) -> Self {
        let mut f = self.clone();
        for m in 0..=self.order() {
            // t -> -t flips the sin terms; the rotation flips y and z.
            f.c[0][m] = self.c[0][m];
            f.s[0][m] = -self.s[0][m];
            f.c[1][m] = -self.c[1][m];
            f.s[1][m] = self.s[1][m];
            f.c[2][m] = -self.c[2][m];
            f.s[2][m] = self.s[2][m];
        }
        f
    }

    /// `n` points at t = 2 pi k / n (the curve is closed; the first point is not repeated).
    pub fn sample(&self, n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|k| self.point(2.0 * PI * k as f64 / n as f64))
            .collect()
    }

    /// Length by the trapezoid rule on |x'(t)| (spectrally accurate for a periodic curve).
    pub fn length(&self, n: usize) -> f64 {
        (0..n)
            .map(|k| {
                let d = self.tangent(2.0 * PI * k as f64 / n as f64);
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .sum::<f64>()
            * 2.0 * PI
            / n as f64
    }
}

/// A toy modular stellarator: `n_coils` planar elliptical coils in poloidal planes, the
/// ellipse turning by pi per field period (a rotating ellipse, l = 2).
#[derive(Clone, Copy, Debug)]
pub struct RotatingEllipse {
    pub n_fp: usize,
    pub coils_per_period: usize,
    /// (m) coil centre major radius.
    pub rc: f64,
    /// (m) ellipse semi-axes.
    pub a_long: f64,
    pub a_short: f64,
    /// (m) coil centres follow a helix about R = rc: (R, Z) = (rc + shift cos(n_fp phi),
    /// shift_z sin(n_fp phi)); shift_z = -shift reverses the helix. Zero gives a planar axis.
    pub shift: f64,
    pub shift_z: f64,
    /// (A-turns) per coil.
    pub current: f64,
}

impl RotatingEllipse {
    pub fn n_coils(&self) -> usize {
        self.n_fp * self.coils_per_period
    }

    /// Toroidal angle of coil k: half a spacing off phi = 0, so phi = 0 and phi = pi / n_fp
    /// fall between coils and are planes of stellarator symmetry.
    pub fn phi(&self, k: usize) -> f64 {
        2.0 * PI * (k as f64 + 0.5) / self.n_coils() as f64
    }

    /// Current per coil for a toroidally averaged field b0 at major radius r0.
    pub fn with_field(mut self, b0: f64, r0: f64) -> Self {
        self.current = 2.0 * PI * r0 * b0 / (MU_0 * self.n_coils() as f64);
        self
    }

    /// Coil k as an order-1 Fourier curve. In the poloidal plane at phi, with
    /// e_R = (cos phi, sin phi, 0) and e_Z = (0, 0, 1), the long axis points along
    /// u = cos(alpha) e_R + sin(alpha) e_Z with alpha = n_fp phi / 2, and
    /// x(t) = c + a_long cos(t) u - a_short sin(t) v, v = -sin(alpha) e_R + cos(alpha) e_Z,
    /// c = (rc + shift cos(n_fp phi)) e_R + shift_z sin(n_fp phi) e_Z.
    /// The minus sign makes positive current give +B_phi inside (u x v = -e_phi).
    pub fn coil(&self, k: usize) -> FourierCurve {
        let phi = self.phi(k);
        let alpha = 0.5 * self.n_fp as f64 * phi;
        let (sp, cp) = phi.sin_cos();
        let (sa, ca) = alpha.sin_cos();
        let er = [cp, sp, 0.0];
        let ez = [0.0, 0.0, 1.0];
        let u: [f64; 3] = std::array::from_fn(|i| ca * er[i] + sa * ez[i]);
        let v: [f64; 3] = std::array::from_fn(|i| -sa * er[i] + ca * ez[i]);
        let (sn, cn) = (self.n_fp as f64 * phi).sin_cos();
        let mut f = FourierCurve::zeros(1);
        for i in 0..3 {
            f.c[i][0] = (self.rc + self.shift * cn) * er[i] + self.shift_z * sn * ez[i];
            f.c[i][1] = self.a_long * u[i];
            f.s[i][1] = -self.a_short * v[i];
        }
        f
    }

    pub fn coils(&self) -> Vec<FourierCurve> {
        (0..self.n_coils()).map(|k| self.coil(k)).collect()
    }
}

/// Coils discretized into straight filaments, laid out for cfsem (structure of arrays).
#[derive(Clone, Debug, Default)]
pub struct Filaments {
    pub x: [Vec<f64>; 3],
    pub dl: [Vec<f64>; 3],
    pub i: Vec<f64>,
    pub r: Vec<f64>,
}

impl Filaments {
    /// Each coil as a closed polygon of `segs` segments through points on the curve.
    pub fn new(coils: &[FourierCurve], current: f64, segs: usize, wire_radius: f64) -> Self {
        let mut f = Filaments::default();
        for c in coils {
            let pts = c.sample(segs);
            for k in 0..segs {
                let a = pts[k];
                let b = pts[(k + 1) % segs];
                for d in 0..3 {
                    f.x[d].push(a[d]);
                    f.dl[d].push(b[d] - a[d]);
                }
                f.i.push(current);
                f.r.push(wire_radius);
            }
        }
        f
    }

    pub fn len(&self) -> usize {
        self.i.len()
    }

    /// B (T) at one Cartesian point.
    pub fn b(&self, p: [f64; 3]) -> [f64; 3] {
        let (mut bx, mut by, mut bz) = ([0.0], [0.0], [0.0]);
        flux_density_linear_filament(
            (&[p[0]], &[p[1]], &[p[2]]),
            (&self.x[0], &self.x[1], &self.x[2]),
            (&self.dl[0], &self.dl[1], &self.dl[2]),
            &self.i,
            &self.r,
            (&mut bx, &mut by, &mut bz),
        )
        .expect("filament arrays line up");
        [bx[0], by[0], bz[0]]
    }

    /// (B_R, B_phi, B_Z) at cylindrical (R, phi, Z).
    pub fn b_cyl(&self, r: f64, phi: f64, z: f64) -> [f64; 3] {
        let (s, c) = phi.sin_cos();
        let b = self.b([r * c, r * s, z]);
        [b[0] * c + b[1] * s, -b[0] * s + b[1] * c, b[2]]
    }
}

/// Minimum distance between two point sets (coil-to-coil clearance from samples).
pub fn min_distance(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let mut d2 = f64::INFINITY;
    for p in a {
        for q in b {
            let e = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
            d2 = d2.min(e);
        }
    }
    d2.sqrt()
}
