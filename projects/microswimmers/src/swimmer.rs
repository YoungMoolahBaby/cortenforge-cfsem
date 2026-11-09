//! The swimmer: a printed helix with a head on its axis holding a magnet magnetized across the
//! axis. A thin strut joins them (its drag ignored). The head is a capsule (`swim`'s desk
//! swimmer) or a "bullet" with a flat, open nose (`shape`: the magnet goes in through the nose).

use std::f64::consts::PI;
use std::time::Instant;

use cfsem::MU_0;
use nalgebra::Matrix6;

use crate::helix::{Helix, V3, bullet_points, capsule_points};
use crate::stokeslet::{self, EPS_C, extrapolate_n};

/// CG tolerance (relative residual).
const TOL: f64 = 1e-10;

#[derive(Clone, Copy, Debug)]
pub struct Swimmer {
    pub helix: Helix,
    /// (m) the head's radius...
    pub head_radius: f64,
    /// (m) ...and straight length; below it a hemisphere, above it a hemisphere or a flat end.
    pub head_len: f64,
    pub flat_nose: bool,
    /// (m) between the helix's end cap and the head.
    pub gap: f64,
    /// (m) the magnet: a cylinder magnetized across its axis, at the bottom of the head's
    /// straight part.
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

/// The model's resistance matrices about `Swimmer::centre`: at the two resolutions, and
/// extrapolated (the one to use).
pub struct Model {
    pub r: Matrix6<f64>,
    pub r12: Matrix6<f64>,
    pub r16: Matrix6<f64>,
    pub points: usize,
    pub seconds: f64,
}

impl Swimmer {
    /// The first desk design: a 3 mm helix of 0.6 mm wire (a resin print), 3 turns, and a
    /// 2 x 2 mm magnet (N42-class, Br 1.3 T as in magnet-tracking) in a 2.5 mm capsule head.
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
            flat_nose: false,
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
    pub fn head_z0(&self) -> f64 {
        let h = &self.helix;
        h.c() * h.t_max() + h.wire + self.gap
    }
    /// z of the nose.
    pub fn nose_z(&self) -> f64 {
        self.head_z0()
            + self.head_radius
            + self.head_len
            + if self.flat_nose { 0.0 } else { self.head_radius }
    }
    /// The magnet's centre, at the bottom of the head's straight part.
    pub fn magnet_centre(&self) -> V3 {
        V3::new(
            0.0,
            0.0,
            self.head_z0() + self.head_radius + 0.5 * self.mag_len,
        )
    }
    /// The head's volume and centroid.
    fn head_volume(&self) -> (f64, V3) {
        let (r, l) = (self.head_radius, self.head_len);
        let base = self.head_z0() + r;
        let cyl = PI * r * r * l;
        let hemi = 2.0 / 3.0 * PI * r.powi(3);
        let (mut v, mut m) = (cyl + hemi, cyl * (base + 0.5 * l) + hemi * (base - 3.0 / 8.0 * r));
        if !self.flat_nose {
            v += hemi;
            m += hemi * (base + l + 3.0 / 8.0 * r);
        }
        (v, V3::new(0.0, 0.0, m / v))
    }
    /// (z of the tail, z of the nose).
    fn ends(&self) -> (f64, f64) {
        let h = &self.helix;
        (-h.wire * h.c() / h.speed(), self.nose_z())
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
        let (head, ch) = self.head_volume();
        let mag = PI * self.mag_radius.powi(2) * self.mag_len;
        // Whole turns: the wire's centroid is on the axis, halfway along.
        let (cw, cm) = (h.centre(), self.magnet_centre());
        let volume = wire + head;
        let mass = self.rho_body * (volume - mag) + self.rho_magnet * mag;
        let com = (cw * (self.rho_body * wire) + ch * (self.rho_body * head)
            - cm * (self.rho_body * mag)
            + cm * (self.rho_magnet * mag))
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
        let (r, l, z0) = (self.head_radius, self.head_len, self.head_z0());
        pts.extend(if self.flat_nose {
            bullet_points(r, l, z0, h)
        } else {
            capsule_points(r, l, z0, h)
        });
        (pts, h)
    }

    /// The regularized-Stokeslet model at 12 and 16 round, extrapolated.
    pub fn model(&self, mu: f64) -> Model {
        self.model_at(mu, 12, 16)
    }

    /// The model at `lo` and `hi` points round the wire, extrapolated (fields `r12`, `r16` hold
    /// the two resolutions).
    pub fn model_at(&self, mu: f64, lo: usize, hi: usize) -> Model {
        let t0 = Instant::now();
        let c = self.centre();
        let solve = |n| {
            let (pts, h) = self.points(n);
            let (r, _) = stokeslet::resistance(&pts, EPS_C * h, mu, &c, TOL);
            (r, pts.len())
        };
        let (r_lo, _) = solve(lo);
        let (r_hi, points) = solve(hi);
        Model {
            r: extrapolate_n(&r_lo, lo, &r_hi, hi),
            r12: r_lo,
            r16: r_hi,
            points,
            seconds: t0.elapsed().as_secs_f64(),
        }
    }
}
