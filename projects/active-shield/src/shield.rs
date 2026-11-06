//! A confined toroidal field shield round a torus habitat: two windings of `n` poloidal loops
//! round the tube (inner at minor radius r1 carrying -I, outer at r2 carrying +I). Between them
//! the field is toroidal, B = mu0 N I / (2 pi R); inside r1 (the habitat) and outside r2 it
//! cancels. Particles from outside must cross the field shell, bending power about B (r2 - r1).

use std::f64::consts::PI;

use cfsem::MU_0;
use cfsem::physics::circular_filament::{
    flux_density_circular_filament_cartesian_par, vector_potential_circular_filament_cartesian_par,
};

/// Geometric mean distance of a thin strip of width w from itself is w e^(-3/2): a row of round
/// filaments that stand for a current sheet gets that equivalent radius (Rosa & Grover).
pub const STRIP_GMD: f64 = 0.223_130_160_148_429_83;

#[derive(Clone, Debug)]
pub struct Shield {
    /// Tube centre's distance from the spin axis (m).
    pub r_major: f64,
    /// Inner and outer winding minor radii (m).
    pub r1: f64,
    pub r2: f64,
    /// Loops per winding.
    pub n: usize,
    /// Ampere-turns of each winding, N I (A).
    pub amp_turns: f64,
}

impl Shield {
    /// A shield whose field at the tube centre's radius is `b_centre` (T).
    pub fn with_field(r_major: f64, r1: f64, r2: f64, n: usize, b_centre: f64) -> Self {
        Self { r_major, r1, r2, n, amp_turns: 2.0 * PI * r_major * b_centre / MU_0 }
    }

    pub fn current(&self) -> f64 {
        self.amp_turns / self.n as f64
    }

    /// Field at radius R from the spin axis inside the shell (T).
    pub fn b_at(&self, big_r: f64) -> f64 {
        MU_0 * self.amp_turns / (2.0 * PI * big_r)
    }

    /// Bending power at the tube centre's radius (T m).
    pub fn bending_power(&self) -> f64 {
        self.b_at(self.r_major) * (self.r2 - self.r1)
    }

    /// Distance from the tube's centre circle (the minor radius coordinate).
    pub fn minor(&self, p: [f64; 3]) -> f64 {
        let big_r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        ((big_r - self.r_major).powi(2) + p[2] * p[2]).sqrt()
    }

    /// The field of continuous windings (exact for N -> infinity).
    pub fn b_ideal(&self, p: [f64; 3]) -> [f64; 3] {
        let rho = self.minor(p);
        if rho <= self.r1 || rho >= self.r2 {
            return [0.0; 3];
        }
        let big_r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        let b = self.b_at(big_r);
        [-b * p[1] / big_r, b * p[0] / big_r, 0.0]
    }

    /// Stored energy of the continuous windings: integral of B^2 / 2 mu0 over the shell,
    /// (mu0 (N I)^2 / 2) (sqrt(R^2 - r1^2) - sqrt(R^2 - r2^2)).
    pub fn energy_ideal(&self) -> f64 {
        let r = self.r_major;
        0.5 * MU_0 * self.amp_turns.powi(2)
            * ((r * r - self.r1 * self.r1).sqrt() - (r * r - self.r2 * self.r2).sqrt())
    }

    /// The discrete loops, for cfsem.
    pub fn loops(&self) -> Loops {
        let mut l = Loops::default();
        let spacing = 2.0 * PI * self.r_major / self.n as f64;
        for (radius, current) in [(self.r1, -self.current()), (self.r2, self.current())] {
            for k in 0..self.n {
                let phi = 2.0 * PI * k as f64 / self.n as f64;
                l.push(
                    radius,
                    current,
                    [self.r_major * phi.cos(), self.r_major * phi.sin(), 0.0],
                    [-phi.sin(), phi.cos(), 0.0],
                    STRIP_GMD * spacing,
                );
            }
        }
        l
    }
}

/// Circular loops in structure-of-arrays form for cfsem's oriented circular-filament kernels.
#[derive(Clone, Debug, Default)]
pub struct Loops {
    pub radius: Vec<f64>,
    pub current: Vec<f64>,
    pub c: [Vec<f64>; 3],
    pub normal: [Vec<f64>; 3],
    pub wire: Vec<f64>,
}

impl Loops {
    pub fn push(&mut self, radius: f64, current: f64, c: [f64; 3], normal: [f64; 3], wire: f64) {
        self.radius.push(radius);
        self.current.push(current);
        for i in 0..3 {
            self.c[i].push(c[i]);
            self.normal[i].push(normal[i]);
        }
        self.wire.push(wire);
    }

    fn eval(&self, pts: &[[f64; 3]], potential: bool) -> Vec<[f64; 3]> {
        let (x, y, z): (Vec<f64>, Vec<f64>, Vec<f64>) =
            (pts.iter().map(|p| p[0]).collect(), pts.iter().map(|p| p[1]).collect(), pts.iter().map(|p| p[2]).collect());
        let (mut bx, mut by, mut bz) = (vec![0.0; pts.len()], vec![0.0; pts.len()], vec![0.0; pts.len()]);
        let args = (
            (&self.radius[..], &self.current[..]),
            (&self.c[0][..], &self.c[1][..], &self.c[2][..]),
            (&self.normal[0][..], &self.normal[1][..], &self.normal[2][..]),
            &self.wire[..],
            (&x[..], &y[..], &z[..]),
        );
        let out = (&mut bx[..], &mut by[..], &mut bz[..]);
        if potential {
            vector_potential_circular_filament_cartesian_par(args.0, args.1, args.2, args.3, args.4, out)
        } else {
            flux_density_circular_filament_cartesian_par(args.0, args.1, args.2, args.3, args.4, out)
        }
        .expect("cfsem circular filament");
        (0..pts.len()).map(|i| [bx[i], by[i], bz[i]]).collect()
    }

    /// Flux density (T) at the points.
    pub fn b(&self, pts: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.eval(pts, false)
    }

    /// Vector potential (Wb/m) at the points.
    pub fn a(&self, pts: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.eval(pts, true)
    }
}

/// A point on the loop of minor radius `r` at toroidal angle 0 and poloidal angle `theta` (from
/// the outer equator), and the loop's tangent there for current circulating right-handed about
/// the loop normal (+y at phi = 0).
pub fn loop_point(r_major: f64, r: f64, theta: f64) -> ([f64; 3], [f64; 3]) {
    let p = [r_major + r * theta.cos(), 0.0, r * theta.sin()];
    let t = [theta.sin(), 0.0, -theta.cos()];
    (p, t)
}

pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Flux through the loop of minor radius `r` at phi = 0 (Wb): the line integral of A round it.
pub fn loop_flux(loops: &Loops, r_major: f64, r: f64, samples: usize) -> f64 {
    let pts: Vec<([f64; 3], [f64; 3])> =
        (0..samples).map(|j| loop_point(r_major, r, 2.0 * PI * (j as f64 + 0.5) / samples as f64)).collect();
    let a = loops.a(&pts.iter().map(|p| p.0).collect::<Vec<_>>());
    let ds = 2.0 * PI * r / samples as f64;
    pts.iter().zip(&a).map(|((_, t), a)| dot(*a, *t) * ds).sum()
}
