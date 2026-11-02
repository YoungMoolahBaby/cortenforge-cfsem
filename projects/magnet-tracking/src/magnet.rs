//! Forward models of an axially magnetized cylindrical magnet: a point dipole (cfsem), current
//! loops on its side (cfsem), and Derby and Olbert's closed form (independent).

use std::f64::consts::PI;

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_density_circular_filament_cartesian;
use cfsem::physics::point_source::dipole::flux_density_dipole_scalar;
use nalgebra::Vector3;

pub type V3 = Vector3<f64>;

/// A cylinder magnet: radius, length (m), remanence Br (T).
#[derive(Clone, Copy, Debug)]
pub struct Magnet {
    pub radius: f64,
    pub length: f64,
    pub br: f64,
}

impl Magnet {
    /// (A m^2) dipole moment: M V with M = Br / mu0.
    pub fn moment(&self) -> f64 {
        self.br / MU_0 * PI * self.radius * self.radius * self.length
    }

    /// (T) field at `obs` of the magnet centred at `pos` with axis `u`, as a cfsem point dipole.
    pub fn dipole(&self, pos: &V3, u: &V3, obs: &V3) -> V3 {
        let m = u * self.moment();
        let b = flux_density_dipole_scalar(
            (pos.x, pos.y, pos.z),
            (m.x, m.y, m.z),
            0.0,
            (obs.x, obs.y, obs.z),
        );
        V3::new(b.0, b.1, b.2)
    }

    /// (T) fields at many points from `n` cfsem current loops on the magnet's side (each Br L /
    /// (mu0 n) of the surface current M = Br / mu0, at the centres of n equal slices).
    pub fn loops(&self, pos: &V3, u: &V3, n: usize, obs: &[V3]) -> Vec<V3> {
        let i = self.br * self.length / (MU_0 * n as f64);
        let centres: Vec<V3> = (0..n)
            .map(|k| pos + u * (-0.5 * self.length + (k as f64 + 0.5) * self.length / n as f64))
            .collect();
        let col = |f: fn(&V3) -> f64, v: &[V3]| -> Vec<f64> { v.iter().map(f).collect() };
        let (r, cur) = (vec![self.radius; n], vec![i; n]);
        let (lx, ly, lz) = (
            col(|v| v.x, &centres),
            col(|v| v.y, &centres),
            col(|v| v.z, &centres),
        );
        let (nx, ny, nz) = (vec![u.x; n], vec![u.y; n], vec![u.z; n]);
        let (ox, oy, oz) = (col(|v| v.x, obs), col(|v| v.y, obs), col(|v| v.z, obs));
        let (mut bx, mut by, mut bz) = (
            vec![0.0; obs.len()],
            vec![0.0; obs.len()],
            vec![0.0; obs.len()],
        );
        flux_density_circular_filament_cartesian(
            (&r, &cur),
            (&lx, &ly, &lz),
            (&nx, &ny, &nz),
            &vec![0.0; n],
            (&ox, &oy, &oz),
            (&mut bx, &mut by, &mut bz),
        )
        .expect("lengths");
        (0..obs.len())
            .map(|k| V3::new(bx[k], by[k], bz[k]))
            .collect()
    }

    /// (T) Derby and Olbert's closed form (Am. J. Phys. 78, 229 (2010)) for a finite solenoid of
    /// surface current Br / mu0, in the magnet's frame: (B_rho, B_z) at radius rho and height z from
    /// the centre.
    pub fn derby(&self, rho: f64, z: f64) -> (f64, f64) {
        let (a, b) = (self.radius, 0.5 * self.length);
        let b0 = self.br / PI;
        let (zp, zm) = (z + b, z - b);
        let alpha = |zz: f64| a / (zz * zz + (rho + a) * (rho + a)).sqrt();
        let beta = |zz: f64| zz / (zz * zz + (rho + a) * (rho + a)).sqrt();
        let gamma = (a - rho) / (a + rho);
        let k = |zz: f64| {
            ((zz * zz + (a - rho) * (a - rho)) / (zz * zz + (a + rho) * (a + rho))).sqrt()
        };
        let b_rho =
            b0 * (alpha(zp) * cel(k(zp), 1.0, 1.0, -1.0) - alpha(zm) * cel(k(zm), 1.0, 1.0, -1.0));
        let b_z = b0 * a / (a + rho)
            * (beta(zp) * cel(k(zp), gamma * gamma, 1.0, gamma)
                - beta(zm) * cel(k(zm), gamma * gamma, 1.0, gamma));
        (b_rho, b_z)
    }

    /// (T) the closed form at a world point, for a magnet at `pos` with axis `u`.
    pub fn exact(&self, pos: &V3, u: &V3, obs: &V3) -> V3 {
        let d = obs - pos;
        let z = d.dot(u);
        let radial = d - u * z;
        let rho = radial.norm();
        let (br, bz) = self.derby(rho, z);
        let rhat = if rho > 0.0 { radial / rho } else { V3::zeros() };
        rhat * br + u * bz
    }
}

/// Bulirsch's generalized complete elliptic integral cel(kc, p, c, s) (as in Derby and Olbert's
/// appendix), converged to 1e-14.
pub fn cel(kc: f64, p: f64, c: f64, s: f64) -> f64 {
    assert!(kc != 0.0, "cel: kc = 0");
    let errtol = 1e-14;
    let mut k = kc.abs();
    let (mut pp, mut cc, mut ss) = (p, c, s);
    let mut em = 1.0;
    if p > 0.0 {
        pp = p.sqrt();
        ss = s / pp;
    } else {
        let mut f = kc * kc;
        let mut q = 1.0 - f;
        let g = 1.0 - pp;
        f -= pp;
        q *= ss - c * pp;
        pp = (f / g).sqrt();
        cc = (c - ss) / g;
        ss = -q / (g * g * pp) + cc * pp;
    }
    let mut f = cc;
    cc += ss / pp;
    let mut g = k / pp;
    ss = 2.0 * (ss + f * g);
    pp += g;
    g = em;
    em += k;
    let mut kk = k;
    while (g - k).abs() > g * errtol {
        k = 2.0 * kk.sqrt();
        kk = k * em;
        f = cc;
        cc += ss / pp;
        g = kk / pp;
        ss = 2.0 * (ss + f * g);
        pp += g;
        g = em;
        em += k;
    }
    0.5 * PI * (ss + cc * em) / (em * (em + pp))
}
