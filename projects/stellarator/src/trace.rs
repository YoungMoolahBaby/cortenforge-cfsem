//! Field-line tracing with the toroidal angle as the time variable:
//! dR/dphi = R B_R / B_phi, dZ/dphi = R B_Z / B_phi, by RK4.
//!
//! The magnetic axis is the fixed point of the one-period map. A line and the axis are traced
//! together, so the line's poloidal angle about the axis is known at every step; the rotational
//! transform is the weighted Birkhoff average of the per-period angle advance (super-convergent
//! on a magnetic surface, Das et al. 2017), with the plain average shown beside it.

use std::f64::consts::PI;

use crate::coils::Filaments;

#[derive(Clone, Copy, Debug)]
pub struct Tracer {
    pub n_fp: usize,
    /// RK4 steps per field period (even, so the half-period plane is hit exactly).
    pub steps: usize,
    /// (m) centre of the coil set in (R, Z) and the radius beyond which a line is lost.
    pub centre: (f64, f64),
    pub lost_radius: f64,
}

pub type Rz = (f64, f64);

impl Tracer {
    pub fn period(&self) -> f64 {
        2.0 * PI / self.n_fp as f64
    }

    fn rhs(f: &Filaments, phi: f64, y: Rz) -> Rz {
        let b = f.b_cyl(y.0, phi, y.1);
        (y.0 * b[0] / b[1], y.0 * b[2] / b[1])
    }

    pub fn rk4(f: &Filaments, phi: f64, y: Rz, h: f64) -> Rz {
        let k1 = Self::rhs(f, phi, y);
        let k2 = Self::rhs(f, phi + 0.5 * h, (y.0 + 0.5 * h * k1.0, y.1 + 0.5 * h * k1.1));
        let k3 = Self::rhs(f, phi + 0.5 * h, (y.0 + 0.5 * h * k2.0, y.1 + 0.5 * h * k2.1));
        let k4 = Self::rhs(f, phi + h, (y.0 + h * k3.0, y.1 + h * k3.1));
        (
            y.0 + h / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0),
            y.1 + h / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1),
        )
    }

    /// One field period from phi = 0.
    pub fn period_map(&self, f: &Filaments, mut y: Rz) -> Rz {
        let h = self.period() / self.steps as f64;
        for k in 0..self.steps {
            y = Self::rk4(f, k as f64 * h, y, h);
        }
        y
    }

    /// Newton on the period map: the magnetic axis where it crosses phi = 0.
    pub fn find_axis(&self, f: &Filaments, guess: Rz) -> (Rz, f64) {
        let mut y = guess;
        let eps = 1e-6;
        let mut res = f64::INFINITY;
        for _ in 0..20 {
            let m = self.period_map(f, y);
            let g = (m.0 - y.0, m.1 - y.1);
            res = g.0.hypot(g.1);
            if res < 1e-13 {
                break;
            }
            let mr = self.period_map(f, (y.0 + eps, y.1));
            let mz = self.period_map(f, (y.0, y.1 + eps));
            // Jacobian of g = M(y) - y.
            let j = [
                [(mr.0 - m.0) / eps - 1.0, (mz.0 - m.0) / eps],
                [(mr.1 - m.1) / eps, (mz.1 - m.1) / eps - 1.0],
            ];
            let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
            let d = (
                (j[1][1] * g.0 - j[0][1] * g.1) / det,
                (-j[1][0] * g.0 + j[0][0] * g.1) / det,
            );
            y = (y.0 - d.0, y.1 - d.1);
        }
        (y, res)
    }
    /// The period map's Jacobian at `y` by central differences.
    pub fn jacobian(&self, f: &Filaments, y: Rz) -> [[f64; 2]; 2] {
        let eps = 1e-6;
        let rp = self.period_map(f, (y.0 + eps, y.1));
        let rm = self.period_map(f, (y.0 - eps, y.1));
        let zp = self.period_map(f, (y.0, y.1 + eps));
        let zm = self.period_map(f, (y.0, y.1 - eps));
        [
            [(rp.0 - rm.0) / (2.0 * eps), (zp.0 - zm.0) / (2.0 * eps)],
            [(rp.1 - rm.1) / (2.0 * eps), (zp.1 - zm.1) / (2.0 * eps)],
        ]
    }

    /// Fixed points of the period map from several starting guesses about `centre`, each
    /// classified by its Jacobian, nearest `centre` first.
    pub fn fixed_points(&self, f: &Filaments, centre: Rz) -> Vec<FixedPoint> {
        let mut out: Vec<FixedPoint> = Vec::new();
        for d in [0.0, -0.02, 0.02, -0.04, 0.04, -0.07, 0.07] {
            let (y, res) = self.find_axis(f, (centre.0 + d, centre.1));
            if !(res < 1e-10) || out.iter().any(|p| (p.rz.0 - y.0).hypot(p.rz.1 - y.1) < 1e-6) {
                continue;
            }
            let j = self.jacobian(f, y);
            out.push(FixedPoint {
                rz: y,
                residual: res,
                trace: j[0][0] + j[1][1],
                det: j[0][0] * j[1][1] - j[0][1] * j[1][0],
            });
        }
        let dist = |p: &FixedPoint| (p.rz.0 - centre.0).hypot(p.rz.1 - centre.1);
        out.sort_by(|a, b| dist(a).total_cmp(&dist(b)));
        out
    }

    /// The magnetic axis: the elliptic fixed point nearest `centre`.
    pub fn axis(&self, f: &Filaments, centre: Rz) -> Option<FixedPoint> {
        self.fixed_points(f, centre).into_iter().find(|p| p.elliptic())
    }

    /// The axis through one period: entry k is the axis after k + 1 steps (the last is the axis
    /// on phi = 0 again, the discrete map's fixed point).
    pub fn axis_path(&self, f: &Filaments, axis: Rz) -> Vec<Rz> {
        let h = self.period() / self.steps as f64;
        let mut a = axis;
        (0..self.steps)
            .map(|k| {
                a = Self::rk4(f, k as f64 * h, a, h);
                a
            })
            .collect()
    }

    /// Trace a line from `start` on phi = 0 for `periods` field periods, measuring its poloidal
    /// angle about the axis (`axis_path`) at every step.
    pub fn trace(&self, f: &Filaments, start: Rz, path: &[Rz], periods: usize) -> Line {
        let h = self.period() / self.steps as f64;
        let mut y = start;
        let a = path[self.steps - 1];
        let mut theta = (y.1 - a.1).atan2(y.0 - a.0);
        let mut line = Line {
            plane0: vec![start],
            plane_half: Vec::new(),
            dtheta: Vec::new(),
            lost_after: None,
            rho_min: f64::INFINITY,
            rho_max: 0.0,
        };
        let mut theta_period = theta;
        'outer: for p in 0..periods {
            for k in 0..self.steps {
                let phi = (p * self.steps + k) as f64 * h;
                y = Self::rk4(f, phi, y, h);
                let rho = (y.0 - self.centre.0).hypot(y.1 - self.centre.1);
                if !(rho < self.lost_radius) {
                    line.lost_after = Some(p);
                    break 'outer;
                }
                let a = path[k];
                let th = (y.1 - a.1).atan2(y.0 - a.0);
                let mut d = th - theta.rem_euclid(2.0 * PI);
                d = (d + PI).rem_euclid(2.0 * PI) - PI;
                theta += d;
                let r = (y.0 - a.0).hypot(y.1 - a.1);
                line.rho_min = line.rho_min.min(r);
                line.rho_max = line.rho_max.max(r);
                if k + 1 == self.steps / 2 {
                    line.plane_half.push(y);
                }
            }
            line.plane0.push(y);
            line.dtheta.push(theta - theta_period);
            theta_period = theta;
        }
        line
    }
}

#[derive(Clone, Debug)]
pub struct Line {
    /// Crossings of phi = 0 (mod one period) and of the half-period plane.
    pub plane0: Vec<Rz>,
    pub plane_half: Vec<Rz>,
    /// Poloidal angle advance about the axis in each period.
    pub dtheta: Vec<f64>,
    pub lost_after: Option<usize>,
    /// Distance from the axis, min and max along the line.
    pub rho_min: f64,
    pub rho_max: f64,
}

impl Line {
    /// Rotational transform per toroidal turn: weighted Birkhoff average and plain mean.
    pub fn iota(&self, n_fp: usize) -> (f64, f64) {
        let n = self.dtheta.len();
        if n == 0 {
            return (f64::NAN, f64::NAN);
        }
        let w = |t: f64| {
            if t <= 0.0 || t >= 1.0 {
                0.0
            } else {
                (-1.0 / (t * (1.0 - t))).exp()
            }
        };
        let (mut sw, mut swd) = (0.0, 0.0);
        for (j, d) in self.dtheta.iter().enumerate() {
            let wj = w((j as f64 + 0.5) / n as f64);
            sw += wj;
            swd += wj * d;
        }
        let mean = self.dtheta.iter().sum::<f64>() / n as f64;
        let k = n_fp as f64 / (2.0 * PI);
        (swd / sw * k, mean * k)
    }
}

/// A fixed point of the one-period map.
#[derive(Clone, Copy, Debug)]
pub struct FixedPoint {
    pub rz: Rz,
    pub residual: f64,
    /// Trace and determinant of the period map's Jacobian (det = 1 for a field-line map).
    pub trace: f64,
    pub det: f64,
}

impl FixedPoint {
    pub fn elliptic(&self) -> bool {
        self.trace.abs() < 2.0
    }

    /// Rotational transform on the axis from the Jacobian, iota = n_fp acos(tr / 2) / 2 pi, as a
    /// magnitude (the sign needs the rotation sense) and only if the transform per period is
    /// under 1/2.
    pub fn iota(&self, n_fp: usize) -> f64 {
        n_fp as f64 * (0.5 * self.trace).clamp(-1.0, 1.0).acos() / (2.0 * PI)
    }
}
