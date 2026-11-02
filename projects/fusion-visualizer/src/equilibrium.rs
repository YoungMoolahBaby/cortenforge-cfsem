//! Free-boundary Grad-Shafranov equilibrium on a fixed grid, by Picard iteration.
//!
//! Total flux Psi = Psi_coils + Psi_plasma. Coils are blocks of cfsem circular filaments outside
//! the grid. The plasma current density follows the Jeon (2015) / FreeGS profile
//! J = lambda (beta0 R/R0 + (1 - beta0) R0/R) (1 - psiN^a1)^a2 inside the last closed surface,
//! lambda set by the plasma current. Each step solves Delta* Psi_p = -2 pi mu0 R J with Psi_p on
//! the grid edge from the plasma current itself, summed exactly through cfsem's filament flux.

use std::f64::consts::PI;

use cfsem::MU_0;
use cfsem::physics::circular_filament::{
    flux_circular_filament_par, flux_density_circular_filament_par,
};

use crate::contour::{contours, inside, perimeter};
use crate::gs::{Grid, Solver, field_from_flux, rhs};

/// A rectangular coil block, `n x n` filaments sharing `current` (A-turns, + along phi).
#[derive(Clone, Debug)]
pub struct Coil {
    pub name: &'static str,
    pub r: f64,
    pub z: f64,
    pub w: f64,
    pub h: f64,
    pub n: usize,
    pub current: f64,
}

impl Coil {
    pub fn filaments(&self) -> impl Iterator<Item = (f64, f64, f64)> + '_ {
        let n = self.n;
        (0..n).flat_map(move |a| {
            (0..n).map(move |b| {
                let r = self.r + self.w * ((a as f64 + 0.5) / n as f64 - 0.5);
                let z = self.z + self.h * ((b as f64 + 0.5) / n as f64 - 0.5);
                (r, z, self.current / (n * n) as f64)
            })
        })
    }
}

fn unzip3(f: impl Iterator<Item = (f64, f64, f64)>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mut a, mut b, mut c) = (Vec::new(), Vec::new(), Vec::new());
    for (x, y, z) in f {
        a.push(x);
        b.push(y);
        c.push(z);
    }
    (a, b, c)
}

/// Flux (Wb) of ideal filaments `(r, z, I)` at points, through cfsem.
pub fn filament_flux(src: &[(f64, f64, f64)], r: &[f64], z: &[f64]) -> Vec<f64> {
    let (rf, zf, i) = unzip3(src.iter().copied());
    let mut out = vec![0.0; r.len()];
    flux_circular_filament_par((&rf, &zf, &i), &vec![0.0; rf.len()], (r, z), &mut out)
        .expect("lengths");
    out
}

/// (Br, Bz) of ideal filaments at points, through cfsem.
pub fn filament_field(src: &[(f64, f64, f64)], r: &[f64], z: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let (rf, zf, i) = unzip3(src.iter().copied());
    let (mut br, mut bz) = (vec![0.0; r.len()], vec![0.0; r.len()]);
    flux_density_circular_filament_par(
        (&rf, &zf, &i),
        &vec![0.0; rf.len()],
        (r, z),
        (&mut br, &mut bz),
    )
    .expect("lengths");
    (br, bz)
}

pub fn coil_filaments(coils: &[Coil]) -> Vec<(f64, f64, f64)> {
    coils
        .iter()
        .flat_map(|c| c.filaments().collect::<Vec<_>>())
        .collect()
}

/// Current profile parameters.
#[derive(Clone, Copy, Debug)]
pub struct Profile {
    /// (A) plasma current.
    pub ip: f64,
    /// Share of the current carried by pressure (beta0 in J's R-weighting).
    pub beta0: f64,
    pub a1: f64,
    pub a2: f64,
    /// (m) reference major radius in the profile.
    pub r0: f64,
}

impl Profile {
    fn shape(&self, psin: f64) -> f64 {
        (1.0 - psin.clamp(0.0, 1.0).powf(self.a1)).powf(self.a2)
    }
}

/// What sets the last closed flux surface.
#[derive(Clone, Copy, Debug)]
pub enum Boundary {
    Limiter(f64, f64),
    XPoint(f64, f64),
}

pub struct Equilibrium {
    pub psi: Vec<f64>,
    pub psi_plasma: Vec<f64>,
    pub jphi: Vec<f64>,
    /// Plasma nodes (inside the last closed surface, connected to the axis).
    pub mask: Vec<bool>,
    pub axis: (f64, f64, f64),
    pub psi_b: f64,
    pub boundary: Boundary,
    /// (A/m^2) profile scale.
    pub lambda: f64,
    pub iterations: usize,
    pub converged: bool,
    /// Last relative change of Psi_plasma.
    pub change: f64,
    /// The controlled coil group's final current, if any.
    pub control_current: Option<f64>,
}

/// Bicubic (Catmull-Rom) interpolation of a grid field.
pub fn interp(grid: &Grid, f: &[f64], r: f64, z: f64) -> f64 {
    let fi = (r - grid.rs[0]) / grid.dr();
    let fj = (z - grid.zs[0]) / grid.dz();
    let i = (fi.floor() as isize).clamp(1, grid.nr() as isize - 3);
    let j = (fj.floor() as isize).clamp(1, grid.nz() as isize - 3);
    let (tx, ty) = (fi - i as f64, fj - j as f64);
    let w = |t: f64| {
        [
            0.5 * (-t * t * t + 2.0 * t * t - t),
            0.5 * (3.0 * t * t * t - 5.0 * t * t + 2.0),
            0.5 * (-3.0 * t * t * t + 4.0 * t * t + t),
            0.5 * (t * t * t - t * t),
        ]
    };
    let (wx, wy) = (w(tx), w(ty));
    let mut s = 0.0;
    for (a, wa) in wx.iter().enumerate() {
        for (b, wb) in wy.iter().enumerate() {
            let ii = (i - 1 + a as isize) as usize;
            let jj = (j - 1 + b as isize) as usize;
            s += wa * wb * f[grid.idx(ii, jj)];
        }
    }
    s
}

/// Refine a grid extremum or saddle at node (i, j) by Newton steps on the bicubic interpolant.
fn refine(grid: &Grid, f: &[f64], r: f64, z: f64) -> (f64, f64, f64) {
    let (mut r, mut z) = (r, z);
    let (hr, hz) = (0.25 * grid.dr(), 0.25 * grid.dz());
    for _ in 0..20 {
        let p = |a: f64, b: f64| interp(grid, f, a, b);
        let fr = (p(r + hr, z) - p(r - hr, z)) / (2.0 * hr);
        let fz = (p(r, z + hz) - p(r, z - hz)) / (2.0 * hz);
        let frr = (p(r + hr, z) - 2.0 * p(r, z) + p(r - hr, z)) / (hr * hr);
        let fzz = (p(r, z + hz) - 2.0 * p(r, z) + p(r, z - hz)) / (hz * hz);
        let frz = (p(r + hr, z + hz) - p(r + hr, z - hz) - p(r - hr, z + hz) + p(r - hr, z - hz))
            / (4.0 * hr * hz);
        let det = frr * fzz - frz * frz;
        if det.abs() < 1e-300 {
            break;
        }
        let dr = (fzz * fr - frz * fz) / det;
        let dz = (frr * fz - frz * fr) / det;
        r -= dr.clamp(-grid.dr(), grid.dr());
        z -= dz.clamp(-grid.dz(), grid.dz());
        if dr.abs() < 1e-9 && dz.abs() < 1e-9 {
            break;
        }
    }
    (r, z, interp(grid, f, r, z))
}

/// Saddle points of `f` among nodes inside `region`.
fn saddles(grid: &Grid, f: &[f64], region: &[(f64, f64)]) -> Vec<(f64, f64, f64)> {
    let mut out: Vec<(f64, f64, f64)> = Vec::new();
    let g2 = |i: usize, j: usize| {
        let fr = (f[grid.idx(i + 1, j)] - f[grid.idx(i - 1, j)]) / (2.0 * grid.dr());
        let fz = (f[grid.idx(i, j + 1)] - f[grid.idx(i, j - 1)]) / (2.0 * grid.dz());
        fr * fr + fz * fz
    };
    for i in 2..grid.nr() - 2 {
        for j in 2..grid.nz() - 2 {
            if !inside(region, (grid.rs[i], grid.zs[j])) {
                continue;
            }
            let c = g2(i, j);
            let local_min = (-1i32..=1)
                .flat_map(|a| (-1i32..=1).map(move |b| (a, b)))
                .filter(|&(a, b)| (a, b) != (0, 0))
                .all(|(a, b)| g2((i as i32 + a) as usize, (j as i32 + b) as usize) >= c);
            if !local_min {
                continue;
            }
            let (r, z, v) = refine(grid, f, grid.rs[i], grid.zs[j]);
            // A saddle: the Hessian is indefinite there.
            let (hr, hz) = (grid.dr(), grid.dz());
            let p = |a: f64, b: f64| interp(grid, f, a, b);
            let frr = (p(r + hr, z) - 2.0 * v + p(r - hr, z)) / (hr * hr);
            let fzz = (p(r, z + hz) - 2.0 * v + p(r, z - hz)) / (hz * hz);
            if frr * fzz < 0.0
                && inside(region, (r, z))
                && !out.iter().any(|o| (o.0 - r).hypot(o.1 - z) < 2.0 * hr)
            {
                out.push((r, z, v));
            }
        }
    }
    out
}

/// A machine: coils outside the grid, a vessel outline the plasma must stay inside, and the
/// limiter points that can touch it.
pub struct Machine {
    pub coils: Vec<Coil>,
    pub region: Vec<(f64, f64)>,
    pub limiter: Vec<(f64, f64)>,
    pub control: Option<Control>,
}

/// Feedback on one coil group (all at the same current): each iteration its current is set so
/// the total flux is equal at two points (an isoflux constraint). Exact, as flux is linear in
/// the current. With the vertical-field coils and the inboard and outboard midplane points, it
/// holds the plasma's radial position, as tokamaks do.
#[derive(Clone, Debug)]
pub struct Control {
    pub coils: Vec<usize>,
    pub points: [(f64, f64); 2],
}

pub struct Options {
    pub relax: f64,
    pub tol: f64,
    pub max_iter: usize,
    /// Average J with its mirror image in Z: up-down symmetric equilibria only. The vertical
    /// position is not controlled otherwise, and elongated plasmas are vertically unstable.
    pub symmetric: bool,
    /// Print each iteration.
    pub verbose: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            relax: 0.5,
            tol: 1e-8,
            max_iter: 400,
            symmetric: true,
            verbose: false,
        }
    }
}

pub fn solve(
    grid: &Grid,
    solver: &Solver,
    m: &Machine,
    prof: &Profile,
    opt: &Options,
    start: Option<&Equilibrium>,
) -> Equilibrium {
    let (rn, zn): (Vec<f64>, Vec<f64>) = grid.nodes().map(|n| (n.0, n.1)).unzip();
    // Coil flux: the fixed coils, plus the controlled group per amp.
    let ctrl = m.control.clone();
    let is_ctrl = |k: usize| ctrl.as_ref().is_some_and(|c| c.coils.contains(&k));
    let fixed: Vec<Coil> = m
        .coils
        .iter()
        .enumerate()
        .filter(|(k, _)| !is_ctrl(*k))
        .map(|(_, c)| c.clone())
        .collect();
    let unit: Vec<Coil> = m
        .coils
        .iter()
        .enumerate()
        .filter(|(k, _)| is_ctrl(*k))
        .map(|(_, c)| Coil {
            current: 1.0,
            ..c.clone()
        })
        .collect();
    let psi_fixed = filament_flux(&coil_filaments(&fixed), &rn, &zn);
    let psi_unit = filament_flux(&coil_filaments(&unit), &rn, &zn);
    let mut i_ctrl = ctrl.as_ref().map_or(0.0, |c| m.coils[c.coils[0]].current);
    let mut psi_coil: Vec<f64> = psi_fixed
        .iter()
        .zip(&psi_unit)
        .map(|(a, b)| a + i_ctrl * b)
        .collect();
    let edge: Vec<usize> = (0..grid.len())
        .filter(|&k| grid.is_boundary(k / grid.nz(), k % grid.nz()))
        .collect();
    let (re, ze): (Vec<f64>, Vec<f64>) = edge.iter().map(|&k| (rn[k], zn[k])).unzip();
    let da = grid.dr() * grid.dz();

    // Start: the previous equilibrium, or a parabolic current in the vessel's middle.
    let mut psi_p = match start {
        Some(e) => e.psi_plasma.clone(),
        None => {
            let (rc, zc) = centroid(&m.region);
            let a = 0.5
                * (m.region.iter().map(|p| p.0).fold(f64::MIN, f64::max)
                    - m.region.iter().map(|p| p.0).fold(f64::MAX, f64::min));
            let j0: Vec<f64> = grid.sample(|r, z| {
                let s = ((r - rc).hypot(z - zc) / (0.6 * a)).powi(2);
                if s < 1.0 { 1.0 - s } else { 0.0 }
            });
            let scale = prof.ip / (j0.iter().sum::<f64>() * da);
            let j: Vec<f64> = j0.iter().map(|x| x * scale).collect();
            plasma_flux(grid, solver, &j, &edge, &re, &ze, &rn, &zn)
        }
    };

    let mut eq = None;
    let mut change = f64::INFINITY;
    let mut near = match start {
        Some(e) => (e.axis.0, e.axis.1),
        None => centroid(&m.region),
    };
    for it in 1..=opt.max_iter {
        if let Some(c) = &ctrl {
            let rest: Vec<f64> = psi_fixed.iter().zip(&psi_p).map(|(a, b)| a + b).collect();
            let [p, q] = c.points;
            let d_rest = interp(grid, &rest, p.0, p.1) - interp(grid, &rest, q.0, q.1);
            let d_unit = interp(grid, &psi_unit, p.0, p.1) - interp(grid, &psi_unit, q.0, q.1);
            i_ctrl = -d_rest / d_unit;
            psi_coil = psi_fixed
                .iter()
                .zip(&psi_unit)
                .map(|(a, b)| a + i_ctrl * b)
                .collect();
        }
        let psi: Vec<f64> = psi_coil.iter().zip(&psi_p).map(|(a, b)| a + b).collect();
        // Magnetic axis: a local flux maximum inside the vessel (positive Ip), the one nearest
        // the last axis. Not the global maximum: the coils' flux can peak at the vessel edge.
        let k0 = local_maxima(grid, &psi, &m.region)
            .into_iter()
            .min_by(|&a, &b| {
                let d = |k: usize| (rn[k] - near.0).hypot(zn[k] - near.1);
                d(a).total_cmp(&d(b))
            })
            .expect("no flux maximum inside the vessel");
        let axis = refine(grid, &psi, rn[k0], zn[k0]);
        // Boundary: the limiter point with the highest flux, unless an X-point inside the
        // vessel sits between it and the axis.
        let (mut psi_b, mut boundary) = m
            .limiter
            .iter()
            .map(|&(r, z)| (interp(grid, &psi, r, z), Boundary::Limiter(r, z)))
            .max_by(|a, b| a.0.total_cmp(&b.0))
            .expect("limiter points");
        for (r, z, v) in saddles(grid, &psi, &m.region) {
            if v > psi_b && v < axis.2 {
                psi_b = v;
                boundary = Boundary::XPoint(r, z);
            }
        }
        // Plasma nodes: below the boundary flux, in the vessel, connected to the axis.
        let psin: Vec<f64> = psi
            .iter()
            .map(|p| (axis.2 - p) / (axis.2 - psi_b))
            .collect();
        let mask = flood(grid, &psin, &m.region, k0);
        let mut j: Vec<f64> = (0..grid.len())
            .map(|k| {
                if mask[k] {
                    let r = rn[k];
                    (prof.beta0 * r / prof.r0 + (1.0 - prof.beta0) * prof.r0 / r)
                        * prof.shape(psin[k])
                } else {
                    0.0
                }
            })
            .collect();
        if opt.symmetric {
            let nz = grid.nz();
            for i in 0..grid.nr() {
                for jj in 0..nz / 2 {
                    let (a, b) = (grid.idx(i, jj), grid.idx(i, nz - 1 - jj));
                    let avg = 0.5 * (j[a] + j[b]);
                    j[a] = avg;
                    j[b] = avg;
                }
            }
        }
        let lambda = prof.ip / (j.iter().sum::<f64>() * da);
        for x in &mut j {
            *x *= lambda;
        }
        let new = plasma_flux(grid, solver, &j, &edge, &re, &ze, &rn, &zn);
        let span = (axis.2 - psi_b).abs();
        change = new
            .iter()
            .zip(&psi_p)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
            / span;
        for (p, n) in psi_p.iter_mut().zip(&new) {
            *p += opt.relax * (n - *p);
        }
        near = (axis.0, axis.1);
        if opt.verbose {
            println!(
                "    it {it:>3}: axis ({:.4}, {:+.4}) psi {:.5}, psi_b {:.5} ({:?}), {} plasma nodes, change {change:.2e}",
                axis.0,
                axis.1,
                axis.2,
                psi_b,
                boundary,
                mask.iter().filter(|&&b| b).count()
            );
        }
        let converged = change < opt.tol;
        if converged || it == opt.max_iter {
            let psi: Vec<f64> = psi_coil.iter().zip(&psi_p).map(|(a, b)| a + b).collect();
            eq = Some(Equilibrium {
                psi,

                psi_plasma: psi_p.clone(),
                jphi: j,
                mask,
                axis,
                psi_b,
                boundary,
                lambda,
                iterations: it,
                converged,
                change,
                control_current: ctrl.as_ref().map(|_| i_ctrl),
            });
            break;
        }
    }
    let mut e = eq.expect("at least one iteration");
    e.change = change;
    e
}

/// Plasma flux for a current density: grid solve with the edge from cfsem's filament flux.
#[allow(clippy::too_many_arguments)]
fn plasma_flux(
    grid: &Grid,
    solver: &Solver,
    j: &[f64],
    edge: &[usize],
    re: &[f64],
    ze: &[f64],
    rn: &[f64],
    zn: &[f64],
) -> Vec<f64> {
    let da = grid.dr() * grid.dz();
    let src: Vec<(f64, f64, f64)> = (0..grid.len())
        .filter(|&k| j[k] != 0.0)
        .map(|k| (rn[k], zn[k], j[k] * da))
        .collect();
    let on_edge = filament_flux(&src, re, ze);
    let mut bvals = vec![0.0; grid.len()];
    for (n, &k) in edge.iter().enumerate() {
        bvals[k] = on_edge[n];
    }
    solver.solve(&rhs(grid, j, |r, z| {
        // rhs() asks by position; look the edge value up by node.
        let i = ((r - grid.rs[0]) / grid.dr()).round() as usize;
        let jj = ((z - grid.zs[0]) / grid.dz()).round() as usize;
        bvals[grid.idx(i, jj)]
    }))
}

/// Nodes with psiN < 1 inside the region, 4-connected to node `start`.
fn flood(grid: &Grid, psin: &[f64], region: &[(f64, f64)], start: usize) -> Vec<bool> {
    let mut mask = vec![false; grid.len()];
    let ok = |k: usize| {
        psin[k] < 1.0 && inside(region, (grid.rs[k / grid.nz()], grid.zs[k % grid.nz()]))
    };
    let mut stack = vec![start];
    while let Some(k) = stack.pop() {
        if mask[k] || !ok(k) {
            continue;
        }
        mask[k] = true;
        let (i, j) = (k / grid.nz(), k % grid.nz());
        if i > 0 {
            stack.push(grid.idx(i - 1, j));
        }
        if i + 1 < grid.nr() {
            stack.push(grid.idx(i + 1, j));
        }
        if j > 0 {
            stack.push(grid.idx(i, j - 1));
        }
        if j + 1 < grid.nz() {
            stack.push(grid.idx(i, j + 1));
        }
    }
    mask
}

fn centroid(poly: &[(f64, f64)]) -> (f64, f64) {
    let n = poly.len() as f64;
    (
        poly.iter().map(|p| p.0).sum::<f64>() / n,
        poly.iter().map(|p| p.1).sum::<f64>() / n,
    )
}

/// Shape and global quantities of an equilibrium.
#[derive(Debug)]
pub struct Summary {
    pub lcfs: Vec<(f64, f64)>,
    pub r_geo: f64,
    pub a: f64,
    pub kappa: f64,
    pub delta_upper: f64,
    pub delta_lower: f64,
    pub beta_p: f64,
    pub li: f64,
    pub ip: f64,
    /// (m^3) plasma volume, summed over plasma nodes.
    pub volume: f64,
}

pub fn summarize(grid: &Grid, eq: &Equilibrium, prof: &Profile) -> Summary {
    // The last closed surface: the contour at psi_b enclosing the axis.
    let lcfs = contours(grid, &eq.psi, eq.psi_b * (1.0 - 1e-12) + eq.axis.2 * 1e-12)
        .into_iter()
        .filter(|c| c.len() > 8 && inside(c, (eq.axis.0, eq.axis.1)))
        .min_by(|a, b| perimeter(a).total_cmp(&perimeter(b)))
        .unwrap_or_default();
    let (rmin, rmax) = lcfs
        .iter()
        .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.0), a.1.max(p.0)));
    let top = lcfs
        .iter()
        .copied()
        .fold((0.0, f64::MIN), |a, p| if p.1 > a.1 { p } else { a });
    let bot = lcfs
        .iter()
        .copied()
        .fold((0.0, f64::MAX), |a, p| if p.1 < a.1 { p } else { a });
    let r_geo = 0.5 * (rmin + rmax);
    let a = 0.5 * (rmax - rmin);

    // Volume integrals over plasma nodes, dV = 2 pi R dR dZ.
    let (rn, zn): (Vec<f64>, Vec<f64>) = grid.nodes().map(|n| (n.0, n.1)).unzip();
    let da = grid.dr() * grid.dz();
    let (br, bz) = field_from_flux(grid, &eq.psi);
    // Pressure from the pressure part of J: 2 pi R dp/dPsi = lambda beta0 (R/R0) g(psiN), so
    // p(psiN) = lambda beta0 / (2 pi R0) (Psi_axis - Psi_b) int_psiN^1 g.
    let n_tab = 2000;
    let mut tail = vec![0.0; n_tab + 1]; // int_x^1 g, on x = k / n_tab
    for k in (0..n_tab).rev() {
        let (x0, x1) = (k as f64 / n_tab as f64, (k + 1) as f64 / n_tab as f64);
        tail[k] = tail[k + 1] + 0.5 * (prof.shape(x0) + prof.shape(x1)) / n_tab as f64;
    }
    let p_scale = eq.lambda * prof.beta0 / (2.0 * PI * prof.r0) * (eq.axis.2 - eq.psi_b);
    let (mut vol, mut p_int, mut bp2_int, mut ip) = (0.0, 0.0, 0.0, 0.0);
    for k in 0..grid.len() {
        if !eq.mask[k] {
            continue;
        }
        let dv = 2.0 * PI * rn[k] * da;
        let psin = (eq.axis.2 - eq.psi[k]) / (eq.axis.2 - eq.psi_b);
        let x = (psin.clamp(0.0, 1.0) * n_tab as f64) as usize;
        vol += dv;
        p_int += p_scale * tail[x.min(n_tab)] * dv;
        bp2_int += (br[k] * br[k] + bz[k] * bz[k]) * dv;
        ip += eq.jphi[k] * da;
        let _ = zn[k];
    }
    let bpa = MU_0 * ip / perimeter(&lcfs);
    Summary {
        r_geo,
        a,
        kappa: (top.1 - bot.1) / (2.0 * a),
        delta_upper: (r_geo - top.0) / a,
        delta_lower: (r_geo - bot.0) / a,
        beta_p: 2.0 * MU_0 * (p_int / vol) / (bpa * bpa),
        li: (bp2_int / vol) / (bpa * bpa),
        volume: vol,
        ip,

        lcfs,
    }
}

/// Nodes inside `region` whose flux exceeds all 8 neighbours.
fn local_maxima(grid: &Grid, f: &[f64], region: &[(f64, f64)]) -> Vec<usize> {
    let mut out = Vec::new();
    for i in 1..grid.nr() - 1 {
        for j in 1..grid.nz() - 1 {
            let k = grid.idx(i, j);
            if !inside(region, (grid.rs[i], grid.zs[j])) {
                continue;
            }
            let peak = (-1i32..=1)
                .flat_map(|a| (-1i32..=1).map(move |b| (a, b)))
                .filter(|&(a, b)| (a, b) != (0, 0))
                .all(|(a, b)| f[grid.idx((i as i32 + a) as usize, (j as i32 + b) as usize)] < f[k]);
            if peak {
                out.push(k);
            }
        }
    }
    out
}
