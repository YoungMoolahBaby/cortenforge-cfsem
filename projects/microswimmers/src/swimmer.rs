//! The desk swimmer: a printed helix with a capsule head on its axis holding a magnet magnetized
//! across the axis. A thin strut joins them (its drag ignored).

use std::f64::consts::PI;
use std::time::Instant;

use cfsem::MU_0;
use nalgebra::Matrix6;

use crate::helix::{Helix, V3, capsule_points};
use crate::stokeslet::{self, EPS_C, extrapolate};

/// CG tolerance (relative residual).
const TOL: f64 = 1e-10;

#[derive(Clone, Copy, Debug)]
pub struct Swimmer {
    pub helix: Helix,
    /// (m) the head: a capsule of this radius...
    pub head_radius: f64,
    /// (m) ...and this straight length (plus its round ends).
    pub head_len: f64,
    /// (m) between the helix's end cap and the head.
    pub gap: f64,
    /// (m) the magnet: a cylinder magnetized across its axis, inside the head.
    pub mag_radius: f64,
    pub mag_len: f64,
    /// (T) the magnet's remanence.
    pub br: f64,
    /// (kg/m^3) the printed body and the magnet.
    pub rho_body: f64,
    pub rho_magnet: f64,
}

/// Volume, mass, centre of mass and centre of buoyancy.
pub struct Mass {
    pub volume: f64,
    pub mass: f64,
    pub com: V3,
    pub cob: V3,
}

/// The model's resistance matrices about `Swimmer::centre`: at 12 and 16 points round the wire,
/// and extrapolated (the one to use).
pub struct Model {
    pub r: Matrix6<f64>,
    pub r12: Matrix6<f64>,
    pub r16: Matrix6<f64>,
    pub points: usize,
    pub seconds: f64,
}

impl Swimmer {
    /// The first desk design: a 3 mm helix of 0.6 mm wire (a resin print), 3 turns, and a
    /// 2 x 2 mm magnet (N42-class, Br 1.3 T as in magnet-tracking) in a 2.5 mm head.
    pub fn desk() -> Self {
        Swimmer {
            helix: Helix {
                radius: 1.5e-3,
                pitch: 6e-3,
                turns: 3.0,
                wire: 0.3e-3,
                hand: 1.0,
            },
            head_radius: 1.25e-3,
            head_len: 2.0e-3,
            gap: 0.5e-3,
            mag_radius: 1.0e-3,
            mag_len: 2.0e-3,
            br: 1.3,
            rho_body: 1150.0,
            rho_magnet: 7500.0,
        }
    }

    /// The same shape scaled by `s` (lengths only).
    pub fn scaled(&self, s: f64) -> Self {
        let mut o = *self;
        o.helix.radius *= s;
        o.helix.pitch *= s;
        o.helix.wire *= s;
        o.head_radius *= s;
        o.head_len *= s;
        o.gap *= s;
        o.mag_radius *= s;
        o.mag_len *= s;
        o
    }

    /// z of the head's lower pole.
    fn head_z0(&self) -> f64 {
        let h = &self.helix;
        h.c() * h.t_max() + h.wire + self.gap
    }
    pub fn head_centre(&self) -> V3 {
        V3::new(
            0.0,
            0.0,
            self.head_z0() + self.head_radius + 0.5 * self.head_len,
        )
    }
    /// (z of the tail, z of the nose): the helix's start cap to the head's far pole.
    fn ends(&self) -> (f64, f64) {
        let h = &self.helix;
        (
            -h.wire * h.c() / h.speed(),
            self.head_z0() + 2.0 * self.head_radius + self.head_len,
        )
    }
    pub fn length(&self) -> f64 {
        let (a, b) = self.ends();
        b - a
    }
    /// The middle of the swimmer on its axis: the reference point for torques.
    pub fn centre(&self) -> V3 {
        let (a, b) = self.ends();
        V3::new(0.0, 0.0, 0.5 * (a + b))
    }
    /// (A m^2) the magnet's moment, Br V / mu0.
    pub fn moment(&self) -> f64 {
        self.br * PI * self.mag_radius.powi(2) * self.mag_len / MU_0
    }

    pub fn mass(&self) -> Mass {
        let h = &self.helix;
        let wire = PI * h.wire.powi(2) * h.length() + 4.0 / 3.0 * PI * h.wire.powi(3);
        let head = PI * self.head_radius.powi(2) * self.head_len
            + 4.0 / 3.0 * PI * self.head_radius.powi(3);
        let mag = PI * self.mag_radius.powi(2) * self.mag_len;
        // Whole turns: the wire's centroid is on the axis, halfway along.
        let (cw, ch) = (h.centre(), self.head_centre());
        let volume = wire + head;
        let mass = self.rho_body * (volume - mag) + self.rho_magnet * mag;
        let com = (cw * (self.rho_body * wire)
            + ch * (self.rho_body * (head - mag) + self.rho_magnet * mag))
            / mass;
        let cob = (cw * wire + ch * head) / volume;
        Mass {
            volume,
            mass,
            com,
            cob,
        }
    }

    /// Surface points at `n_around` round the wire, the head at the same spacing.
    pub fn points(&self, n_around: usize) -> (Vec<V3>, f64) {
        let (mut pts, h) = self.helix.tube_points(n_around);
        pts.extend(capsule_points(
            self.head_radius,
            self.head_len,
            self.head_z0(),
            h,
        ));
        (pts, h)
    }

    /// The regularized-Stokeslet model at 12 and 16 round, extrapolated.
    pub fn model(&self, mu: f64) -> Model {
        let t0 = Instant::now();
        let c = self.centre();
        let solve = |n| {
            let (pts, h) = self.points(n);
            let (r, _) = stokeslet::resistance(&pts, EPS_C * h, mu, &c, TOL);
            (r, pts.len())
        };
        let (r12, _) = solve(12);
        let (r16, points) = solve(16);
        Model {
            r: extrapolate(&r12, &r16),
            r12,
            r16,
            points,
            seconds: t0.elapsed().as_secs_f64(),
        }
    }
}
