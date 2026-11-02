//! Toroidal field coils: N circular coils in poloidal planes, through cfsem's oriented circular
//! filament (`flux_density_circular_filament_cartesian_scalar`), cross-checked against a polygon
//! of cfsem's linear filaments.

use std::f64::consts::PI;

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_density_circular_filament_cartesian_scalar;
use cfsem::physics::linear_filament::flux_density_linear_filament_scalar;

/// A ring of identical circular TF coils.
#[derive(Clone, Copy, Debug)]
pub struct TfSet {
    pub n: usize,
    /// (m) coil centre major radius and coil radius.
    pub rc: f64,
    pub radius: f64,
    /// (A-turns) per coil.
    pub current: f64,
}

impl TfSet {
    /// Current per coil for a vacuum field b0 at major radius r0: B = mu0 N I / (2 pi R).
    pub fn for_field(n: usize, rc: f64, radius: f64, b0: f64, r0: f64) -> Self {
        TfSet {
            n,
            rc,
            radius,
            current: 2.0 * PI * r0 * b0 / (MU_0 * n as f64),
        }
    }

    /// R B_phi of the ideal axisymmetric limit (T m).
    pub fn f0(&self) -> f64 {
        MU_0 * self.n as f64 * self.current / (2.0 * PI)
    }

    /// (Bx, By, Bz) at a Cartesian point.
    pub fn field(&self, p: (f64, f64, f64)) -> (f64, f64, f64) {
        let mut b = (0.0, 0.0, 0.0);
        for k in 0..self.n {
            let phi = 2.0 * PI * k as f64 / self.n as f64;
            let (s, c) = phi.sin_cos();
            // Normal along +phi: positive current gives +B_phi inside.
            let d = flux_density_circular_filament_cartesian_scalar(
                (self.radius, self.current),
                (self.rc * c, self.rc * s, 0.0),
                (-s, c, 0.0),
                0.0,
                p,
            );
            b = (b.0 + d.0, b.1 + d.1, b.2 + d.2);
        }
        b
    }

    /// (B_R, B_phi, B_Z) at cylindrical (R, phi, Z).
    pub fn field_cyl(&self, r: f64, phi: f64, z: f64) -> (f64, f64, f64) {
        let (s, c) = phi.sin_cos();
        let (bx, by, bz) = self.field((r * c, r * s, z));
        (bx * c + by * s, -bx * s + by * c, bz)
    }
}

/// One coil (the one at phi = 0) as an `m`-segment polygon of linear filaments.
pub fn coil_polygon_field(tf: &TfSet, m: usize, p: (f64, f64, f64)) -> (f64, f64, f64) {
    // The coil lies in the x-z plane centred at (rc, 0, 0) with normal +y. Walking
    // (x, z) = (rc + a cos t, a sin t) with t falling from 0 to -2 pi puts +B_y at the centre
    // (right-hand rule about +y), matching the circular kernel.
    let pt = |t: f64| (tf.rc + tf.radius * t.cos(), 0.0, tf.radius * t.sin());
    let mut b = (0.0, 0.0, 0.0);
    for k in 0..m {
        let t0 = -2.0 * PI * k as f64 / m as f64;
        let t1 = -2.0 * PI * (k + 1) as f64 / m as f64;
        let d = flux_density_linear_filament_scalar((pt(t0), pt(t1), tf.current), 0.0, p);
        b = (b.0 + d.0, b.1 + d.1, b.2 + d.2);
    }
    b
}

/// The phi = 0 coil alone through the circular kernel.
pub fn coil_circular_field(tf: &TfSet, p: (f64, f64, f64)) -> (f64, f64, f64) {
    flux_density_circular_filament_cartesian_scalar(
        (tf.radius, tf.current),
        (tf.rc, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        0.0,
        p,
    )
}
