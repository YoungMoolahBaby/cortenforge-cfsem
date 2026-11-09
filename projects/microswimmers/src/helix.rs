//! A rigid helix on the z axis, and point clouds on body surfaces for the regularized-Stokeslet
//! model (the helix as a tube with round end caps, spheres and spheroids).

use std::f64::consts::PI;

use nalgebra::Vector3;

pub type V3 = Vector3<f64>;

/// Centreline (R cos t, h R sin t, c t) for t in [0, 2 pi turns], c = pitch / 2 pi. `hand` is +1
/// for a right-handed helix and -1 for its mirror image (y reflected).
#[derive(Clone, Copy, Debug)]
pub struct Helix {
    /// (m) helix radius R, axis to the wire's centreline.
    pub radius: f64,
    /// (m) advance along the axis per turn.
    pub pitch: f64,
    pub turns: f64,
    /// (m) radius of the wire itself, b.
    pub wire: f64,
    pub hand: f64,
}

impl Helix {
    pub fn c(&self) -> f64 {
        self.pitch / (2.0 * PI)
    }
    pub fn t_max(&self) -> f64 {
        2.0 * PI * self.turns
    }
    /// ds/dt along the centreline.
    pub fn speed(&self) -> f64 {
        (self.radius * self.radius + self.c() * self.c()).sqrt()
    }
    pub fn length(&self) -> f64 {
        self.speed() * self.t_max()
    }
    /// Pitch angle between the wire and the axis, degrees (45 when the pitch is 2 pi R).
    pub fn pitch_angle_deg(&self) -> f64 {
        (2.0 * PI * self.radius).atan2(self.pitch).to_degrees()
    }
    /// The middle of the axis: the reference point for torques.
    pub fn centre(&self) -> V3 {
        V3::new(0.0, 0.0, 0.5 * self.c() * self.t_max())
    }
    pub fn point(&self, t: f64) -> V3 {
        V3::new(
            self.radius * t.cos(),
            self.hand * self.radius * t.sin(),
            self.c() * t,
        )
    }
    /// Unit tangent, unit normal (towards the axis; any perpendicular if R = 0) and binormal.
    pub fn frame(&self, t: f64) -> (V3, V3, V3) {
        let tg = V3::new(
            -self.radius * t.sin(),
            self.hand * self.radius * t.cos(),
            self.c(),
        ) / self.speed();
        let n = if self.radius > 0.0 {
            V3::new(-t.cos(), -self.hand * t.sin(), 0.0)
        } else {
            V3::new(1.0, 0.0, 0.0)
        };
        (tg, n, tg.cross(&n))
    }

    /// Points on the wire's surface (a tube of radius b) with hemispherical end caps, about
    /// `n_around` points round the wire and the same spacing along it; rings staggered by half a
    /// spacing. Returns the points and the spacing h.
    pub fn tube_points(&self, n_around: usize) -> (Vec<V3>, f64) {
        let b = self.wire;
        let h = 2.0 * PI * b / n_around as f64;
        let n_rings = (self.length() / h).round().max(1.0) as usize;
        let mut pts = Vec::new();
        for k in 0..=n_rings {
            let t = self.t_max() * k as f64 / n_rings as f64;
            let (_, n, bn) = self.frame(t);
            let x = self.point(t);
            for j in 0..n_around {
                let phi = 2.0 * PI * (j as f64 + 0.5 * (k % 2) as f64) / n_around as f64;
                pts.push(x + (n * phi.cos() + bn * phi.sin()) * b);
            }
        }
        for (t, sign) in [(0.0, -1.0), (self.t_max(), 1.0)] {
            let (tg, n, bn) = self.frame(t);
            let out = tg * sign;
            let x = self.point(t);
            let n_cap = (0.5 * PI * b / h).round().max(1.0) as usize;
            for i in 1..=n_cap {
                // Polar angle from the outward tangent: the rim is pi/2, the pole 0.
                let alpha = 0.5 * PI * (1.0 - i as f64 / n_cap as f64);
                let rho = b * alpha.sin();
                let m = (2.0 * PI * rho / h).round().max(1.0) as usize;
                for j in 0..m {
                    let phi = 2.0 * PI * (j as f64 + 0.5 * (i % 2) as f64) / m as f64;
                    pts.push(x + out * (b * alpha.cos()) + (n * phi.cos() + bn * phi.sin()) * rho);
                }
            }
        }
        (pts, h)
    }
}

/// Points on a spheroid centred at the origin, semi-axis `a` along z and `b` across, in rings at
/// equal arc length (spacing about h) along the meridian, each ring turned by the golden angle.
pub fn spheroid_points(a: f64, b: f64, h: f64) -> Vec<V3> {
    // Arc length along the meridian (b sin u, a cos u), u from 0 (north pole) to pi.
    let n = 20_000;
    let mut s = vec![0.0; n + 1];
    for i in 1..=n {
        let u = PI * (i as f64 - 0.5) / n as f64;
        let ds = (b * b * u.cos().powi(2) + a * a * u.sin().powi(2)).sqrt() * PI / n as f64;
        s[i] = s[i - 1] + ds;
    }
    let total = s[n];
    let n_rings = (total / h).round().max(2.0) as usize;
    let golden = PI * (3.0 - 5f64.sqrt());
    let mut pts = Vec::new();
    let mut i = 0;
    for k in 0..=n_rings {
        let target = total * k as f64 / n_rings as f64;
        while i < n && s[i + 1] < target {
            i += 1;
        }
        let frac = if i < n {
            (target - s[i]) / (s[i + 1] - s[i])
        } else {
            1.0
        };
        let u = PI * (i as f64 + frac) / n as f64;
        let (rho, z) = (b * u.sin(), a * u.cos());
        let m = (2.0 * PI * rho / h).round().max(1.0) as usize;
        for j in 0..m {
            let phi = 2.0 * PI * j as f64 / m as f64 + golden * k as f64;
            pts.push(V3::new(rho * phi.cos(), rho * phi.sin(), z));
        }
    }
    pts
}

/// Points on a capsule on the z axis (radius `r`, straight part `len`), its lower pole at `z0`, in
/// rings at equal arc length (spacing about h) along the meridian, each turned by the golden angle.
pub fn capsule_points(r: f64, len: f64, z0: f64, h: f64) -> Vec<V3> {
    let quarter = 0.5 * PI * r;
    let total = 2.0 * quarter + len;
    let n_rings = (total / h).round().max(2.0) as usize;
    let golden = PI * (3.0 - 5f64.sqrt());
    let mut pts = Vec::new();
    for k in 0..=n_rings {
        let s = total * k as f64 / n_rings as f64;
        let (rho, z) = if s < quarter {
            let a = s / r;
            (r * a.sin(), z0 + r - r * a.cos())
        } else if s <= quarter + len {
            (r, z0 + r + (s - quarter))
        } else {
            let a = (s - quarter - len) / r;
            (r * a.cos(), z0 + r + len + r * a.sin())
        };
        let m = (2.0 * PI * rho / h).round().max(1.0) as usize;
        for j in 0..m {
            let phi = 2.0 * PI * j as f64 / m as f64 + golden * k as f64;
            pts.push(V3::new(rho * phi.cos(), rho * phi.sin(), z));
        }
    }
    pts
}
