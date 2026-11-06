//! Charged-particle tracing through the shield and transmission by backtracing.
//!
//! In a static magnetic field the path depends only on rigidity (momentum per charge), and the
//! direction turns at 1 / (B rho) per metre: the Boris rotation here steps in path length and
//! keeps the speed exactly. By Liouville's theorem the cosmic-ray intensity arriving at a point
//! in a direction is the free-space intensity if the reversed path (opposite charge, opposite
//! direction) escapes to infinity, and zero if it comes back to the habitat or stays trapped.
//! So the transmission is the fraction of arrival directions, weighted by their flux through the
//! habitat's surface, whose reversed paths escape.

use std::f64::consts::PI;

use crate::shield::{cross, dot, norm, Shield};

pub const C: f64 = 299_792_458.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fate {
    /// Reached free space: the direction is open to cosmic rays.
    Escaped,
    /// Came back to the habitat: closed.
    Habitat,
    /// Still in the shield after the path-length limit: closed (and rare).
    Trapped,
}

/// Rotate `d` (unit) by the field `b` over a path step `ds` for a particle of signed rigidity
/// `brho` (T m, B rho = p / q; negative for negative charge): the Boris rotation, exact in speed.
pub fn boris_turn(d: [f64; 3], b: [f64; 3], ds: f64, brho: f64) -> [f64; 3] {
    let f = 0.5 * ds / brho;
    let t = [b[0] * f, b[1] * f, b[2] * f];
    let t2 = dot(t, t);
    let s = [2.0 * t[0] / (1.0 + t2), 2.0 * t[1] / (1.0 + t2), 2.0 * t[2] / (1.0 + t2)];
    let dt = cross(d, t);
    let dp = [d[0] + dt[0], d[1] + dt[1], d[2] + dt[2]];
    let ds2 = cross(dp, s);
    [d[0] + ds2[0], d[1] + ds2[1], d[2] + ds2[2]]
}

/// Settings for one trace.
#[derive(Clone, Copy, Debug)]
pub struct Tracer {
    /// Habitat outer radius (minor radius, m).
    pub r_hab: f64,
    /// Largest step in the field (m) and the most turning per step (rad).
    pub ds_max: f64,
    pub turn_max: f64,
    /// Straight-line march step outside the field (m).
    pub march: f64,
    /// Path-length limit (m) before a path counts as trapped.
    pub max_path: f64,
}

impl Default for Tracer {
    fn default() -> Self {
        Self { r_hab: 70.0, ds_max: 0.25, turn_max: 0.05, march: 0.5, max_path: 60_000.0 }
    }
}

impl Tracer {
    /// Follow a particle of signed rigidity `brho` (T m) from `p` along `d` through the ideal
    /// field of `sh` until it escapes, reaches the habitat or runs out of path.
    pub fn follow(&self, sh: &Shield, mut p: [f64; 3], mut d: [f64; 3], brho: f64) -> Fate {
        // The whole torus lies within this distance of the origin: a straight path outside it
        // heading away never comes back.
        let far = sh.r_major + sh.r2;
        let mut path = 0.0;
        while path < self.max_path {
            let rho = sh.minor(p);
            if rho < sh.r1 {
                // Field-free inside: straight on to the habitat or back into the shell. The
                // minor radius is the distance to the tube's centre circle, so the distances
                // to both surfaces bound a safe straight step (sphere tracing).
                loop {
                    let r = sh.minor(p);
                    if r - self.r_hab <= 1e-7 {
                        return Fate::Habitat;
                    }
                    if sh.r1 - r <= 1e-7 {
                        // Nudge across into the shell.
                        p = [p[0] + 1e-4 * d[0], p[1] + 1e-4 * d[1], p[2] + 1e-4 * d[2]];
                        break;
                    }
                    let s = (sh.r1 - r).min(r - self.r_hab).max(self.march * 1e-3);
                    p = [p[0] + s * d[0], p[1] + s * d[1], p[2] + s * d[2]];
                    path += s;
                    if path > self.max_path {
                        return Fate::Trapped;
                    }
                }
            } else if rho > sh.r2 {
                // Field-free outside: straight on until it re-enters the shell (the ring's far
                // side) or leaves the torus's bounding sphere heading away.
                loop {
                    let r = sh.minor(p);
                    if r - sh.r2 <= 1e-7 {
                        p = [p[0] + 1e-4 * d[0], p[1] + 1e-4 * d[1], p[2] + 1e-4 * d[2]];
                        break;
                    }
                    if dot(p, p).sqrt() > far && dot(p, d) >= 0.0 {
                        return Fate::Escaped;
                    }
                    let s = (r - sh.r2).max(self.march * 1e-3);
                    p = [p[0] + s * d[0], p[1] + s * d[1], p[2] + s * d[2]];
                    path += s;
                    if path > 1e3 * self.max_path {
                        return Fate::Escaped;
                    }
                }
            } else {
                // In the field: Boris steps in path length.
                let b = sh.b_ideal(p);
                let bmag = norm(b).max(1e-30);
                let ds = (self.turn_max * brho.abs() / bmag).min(self.ds_max);
                d = boris_turn(d, b, ds, brho);
                p = [p[0] + ds * d[0], p[1] + ds * d[1], p[2] + ds * d[2]];
                path += ds;
            }
        }
        Fate::Trapped
    }
}

/// Cosine-weighted arrival directions over the inward hemisphere at a surface point with outward
/// normal `n` and tangents `e1`, `e2` (a Fibonacci spiral on the unit disk, projected up):
/// returns the reversed (outward-going) directions to backtrace.
pub fn reversed_arrivals(n: [f64; 3], e1: [f64; 3], e2: [f64; 3], count: usize) -> Vec<[f64; 3]> {
    let golden = PI * (3.0 - 5f64.sqrt());
    (0..count)
        .map(|i| {
            let r = ((i as f64 + 0.5) / count as f64).sqrt(); // sin(alpha)
            let beta = golden * i as f64;
            let c = (1.0 - r * r).sqrt(); // cos(alpha)
            let (sb, cb) = beta.sin_cos();
            [
                c * n[0] + r * (cb * e1[0] + sb * e2[0]),
                c * n[1] + r * (cb * e1[1] + sb * e2[1]),
                c * n[2] + r * (cb * e1[2] + sb * e2[2]),
            ]
        })
        .collect()
}

/// Rigidity (V) to B rho (T m).
pub fn brho(rigidity_v: f64) -> f64 {
    rigidity_v / C
}

/// Transmission of an isotropic cosmic-ray flux of rigidity `rigidity_v` (V, positive charge)
/// through the shield to the habitat's surface: per poloidal angle and area-weighted overall.
pub fn transmission(sh: &Shield, tr: &Tracer, rigidity_v: f64, poloidal: usize, directions: usize) -> (f64, Vec<(f64, f64)>) {
    use rayon::prelude::*;
    let back = -brho(rigidity_v); // reversed charge for backtracing
    let rows: Vec<(f64, f64, f64)> = (0..poloidal)
        .into_par_iter()
        .map(|j| {
            let theta = 2.0 * PI * (j as f64 + 0.5) / poloidal as f64;
            let (st, ct) = theta.sin_cos();
            let n = [ct, 0.0, st];
            let e1 = [0.0, 1.0, 0.0];
            let e2 = [-st, 0.0, ct];
            let p = [sh.r_major + tr.r_hab * ct, 0.0, tr.r_hab * st];
            let dirs = reversed_arrivals(n, e1, e2, directions);
            let open = dirs
                .iter()
                .filter(|d| {
                    let q = [p[0] + 1e-3 * d[0], p[1] + 1e-3 * d[1], p[2] + 1e-3 * d[2]];
                    tr.follow(sh, q, **d, back) == Fate::Escaped
                })
                .count();
            let weight = sh.r_major + tr.r_hab * ct; // area element ~ R
            (theta, open as f64 / directions as f64, weight)
        })
        .collect();
    let wsum: f64 = rows.iter().map(|r| r.2).sum();
    let total = rows.iter().map(|r| r.1 * r.2).sum::<f64>() / wsum;
    (total, rows.into_iter().map(|r| (r.0, r.1)).collect())
}
