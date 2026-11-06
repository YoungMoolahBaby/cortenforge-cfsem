//! The coil search's design space and score.
//!
//! Design: the three coils of a half period (phi = 10, 30, 50 deg) are free; the other fifteen
//! follow by stellarator symmetry (`FourierCurve::mirror`) and the field period
//! (`FourierCurve::rotate_z`), so every candidate keeps both symmetries exactly. Each free coil
//! is the baseline coil plus a displacement in its own frame (e_R, e_phi, e_Z at the coil's
//! angle), a Fourier series to order 2: 5 coefficients a direction, 45 numbers in all, each a
//! fraction in [-1, 1] of a per-order scale. e_phi modes bend the coil out of its plane.
//!
//! Score (higher is better), from field lines traced through the candidate:
//! - r_good: lines start 5, 10, ..., 60 mm out from the axis on the phi = 0 midplane and run 150
//!   periods; a line is good if it stays `CLEARANCE` from every coil and the Birkhoff iota of
//!   its two halves agree to `IOTA_TOL` (a magnetic surface, not a chaotic layer). r_good is the
//!   start of the last good line in an unbroken run from the axis (mm).
//! - plus `IOTA_WEIGHT` x the lowest iota on those good lines, capped at `IOTA_CAP`: more
//!   transform makes surfaces stiffer and is what modular coils are for.
//! - minus penalties for a coil bent tighter than `RC_MIN` or two coils closer than `GAP_MIN`
//!   (centrelines), the limits for winding 8 mm conductor bundles on a printed former.

use std::f64::consts::PI;
use std::time::Instant;

use rayon::prelude::*;

use crate::coils::{Filaments, FourierCurve, RotatingEllipse};
use crate::design;
use crate::trace::{FixedPoint, Rz, Tracer};

pub const ORDER: usize = 2;
pub const DIRS: usize = 3;
pub const COEFFS: usize = 2 * ORDER + 1;
pub const UNIQUE: usize = 3;
pub const N_PARAMS: usize = UNIQUE * DIRS * COEFFS;
/// (m) displacement scale by mode order (0, 1, 2).
pub const SCALE: [f64; ORDER + 1] = [0.010, 0.015, 0.008];

/// (m) field lines must stay this far from every coil centreline (4 mm conductor + 6 mm of
/// former and room).
pub const CLEARANCE: f64 = 0.010;
pub const IOTA_TOL: f64 = 1e-3;
pub const IOTA_WEIGHT: f64 = 50.0;
pub const IOTA_CAP: f64 = 0.25;
/// (m) least radius of curvature, and least distance between coil centrelines.
pub const RC_MIN: f64 = 0.012;
pub const GAP_MIN: f64 = 0.030;
/// Field periods traced per line, and line starts (m from the axis).
pub const PERIODS: usize = 300;
pub const LINES: usize = 12;
pub const LINE_STEP: f64 = 0.005;
/// Coarser than `trace` for speed: 64 filaments a coil, 60 RK4 steps a period.
pub const SEGS: usize = 64;
pub const STEPS: usize = 120;

#[derive(Clone, Debug, PartialEq)]
pub struct Design {
    /// N_PARAMS fractions in [-1, 1].
    pub x: Vec<f64>,
}

impl Design {
    pub fn zero() -> Self {
        Design { x: vec![0.0; N_PARAMS] }
    }

    fn idx(j: usize, d: usize, q: usize) -> usize {
        (j * DIRS + d) * COEFFS + q
    }

    /// The three free coils (phi = 10, 30, 50 deg for the baseline).
    pub fn free_coils(&self, base: &RotatingEllipse) -> Vec<FourierCurve> {
        (0..UNIQUE)
            .map(|j| {
                let mut c = base.coil(j).with_order(ORDER);
                let (s, co) = base.phi(j).sin_cos();
                let frame = [[co, s, 0.0], [-s, co, 0.0], [0.0, 0.0, 1.0]];
                for (d, e) in frame.iter().enumerate() {
                    let a0 = self.x[Self::idx(j, d, 0)].clamp(-1.0, 1.0) * SCALE[0];
                    for k in 0..3 {
                        c.c[k][0] += a0 * e[k];
                    }
                    for m in 1..=ORDER {
                        let a = self.x[Self::idx(j, d, 2 * m - 1)].clamp(-1.0, 1.0) * SCALE[m];
                        let b = self.x[Self::idx(j, d, 2 * m)].clamp(-1.0, 1.0) * SCALE[m];
                        for k in 0..3 {
                            c.c[k][m] += a * e[k];
                            c.s[k][m] += b * e[k];
                        }
                    }
                }
                c
            })
            .collect()
    }

    /// All coils: each free coil and its stellarator partner, copied into every field period.
    pub fn coils(&self, base: &RotatingEllipse) -> Vec<FourierCurve> {
        assert_eq!(base.coils_per_period, 2 * UNIQUE, "three free coils a half period");
        let free = self.free_coils(base);
        let per = 2.0 * PI / base.n_fp as f64;
        let mut out = Vec::new();
        for p in 0..base.n_fp {
            for c in &free {
                out.push(c.rotate_z(p as f64 * per));
                out.push(c.mirror().rotate_z(p as f64 * per));
            }
        }
        out
    }
}

#[derive(Clone, Debug, Default)]
pub struct Score {
    pub value: f64,
    /// (mm) start of the outermost good line in the unbroken run from the axis.
    pub r_good: f64,
    pub iota_axis: f64,
    pub iota_min: f64,
    /// (mm) tightest radius of curvature, closest pair of coil centrelines.
    pub rc_min: f64,
    pub gap: f64,
    pub penalty: f64,
    pub axis: Option<Rz>,
    pub note: &'static str,
    pub secs: f64,
}

/// Tightest radius of curvature (m) and closest distance between centrelines of different coils
/// (m), checked for the free coils against all coils (symmetry covers the rest).
pub fn geometry(free: &[FourierCurve], all: &[FourierCurve]) -> (f64, f64) {
    let n = 256;
    let mut kmax = 0.0_f64;
    for c in free {
        for i in 0..n {
            kmax = kmax.max(c.curvature(2.0 * PI * i as f64 / n as f64));
        }
    }
    let samples: Vec<Vec<[f64; 3]>> = all.iter().map(|c| c.sample(128)).collect();
    let mut gap = f64::INFINITY;
    for c in free {
        let s = c.sample(128);
        for other in &samples {
            // Skip the coil itself (its own samples coincide).
            let same = other
                .iter()
                .zip(&s)
                .all(|(p, q)| (p[0] - q[0]).abs() + (p[1] - q[1]).abs() + (p[2] - q[2]).abs() < 1e-12);
            if same {
                continue;
            }
            gap = gap.min(crate::coils::min_distance(&s, other));
        }
    }
    (1.0 / kmax, gap)
}

pub fn tracer(base: &RotatingEllipse) -> Tracer {
    Tracer {
        n_fp: base.n_fp,
        steps: STEPS,
        centre: (base.rc, 0.0),
        lost_radius: base.a_long + base.shift.abs() + 0.03,
    }
}

/// The elliptic fixed point nearest the coil centre. Stellarator symmetry keeps the axis on
/// Z = 0 at phi = 0, so Newton starts from six radii along that line, in parallel.
pub fn find_axis(tr: &Tracer, f: &Filaments, base: &RotatingEllipse) -> Option<FixedPoint> {
    let centre = (base.rc, 0.0);
    let found: Vec<(Rz, f64)> = [0.0, 0.02, -0.02, 0.04, -0.04, 0.06]
        .par_iter()
        .map(|d| tr.find_axis(f, (centre.0 + d, centre.1)))
        .collect();
    let mut fps: Vec<FixedPoint> = Vec::new();
    for (y, res) in found {
        if !(res < 1e-10) || fps.iter().any(|p| (p.rz.0 - y.0).hypot(p.rz.1 - y.1) < 1e-6) {
            continue;
        }
        let j = tr.jacobian(f, y);
        fps.push(FixedPoint {
            rz: y,
            residual: res,
            trace: j[0][0] + j[1][1],
            det: j[0][0] * j[1][1] - j[0][1] * j[1][0],
        });
    }
    let dist = |p: &FixedPoint| (p.rz.0 - centre.0).hypot(p.rz.1 - centre.1);
    fps.into_iter()
        .filter(|p| p.elliptic())
        .min_by(|a, b| dist(a).total_cmp(&dist(b)))
}

/// Per-line detail from `evaluate_lines`.
#[derive(Clone, Debug)]
pub struct LineResult {
    pub start: f64,
    pub lost_after: Option<usize>,
    pub iota: f64,
    pub halves: (f64, f64),
    pub good: bool,
}

pub fn evaluate(d: &Design, base: &RotatingEllipse) -> Score {
    evaluate_lines(d, base).0
}

pub fn evaluate_lines(d: &Design, base: &RotatingEllipse) -> (Score, Vec<LineResult>) {
    let t0 = Instant::now();
    let free = d.free_coils(base);
    let all = d.coils(base);
    let (rc, gap) = geometry(&free, &all);
    let penalty = 2e3 * ((RC_MIN - rc).max(0.0) + (GAP_MIN - gap).max(0.0));
    let mut s = Score {
        rc_min: rc * 1e3,
        gap: gap * 1e3,
        penalty,
        ..Score::default()
    };
    let done = |mut s: Score, value: f64, note: &'static str| {
        s.value = value;
        s.note = note;
        s.secs = t0.elapsed().as_secs_f64();
        s
    };
    // Hopeless geometry: don't trace.
    if rc < 0.4 * RC_MIN || gap < 0.4 * GAP_MIN {
        return (done(s, -50.0 - penalty, "geometry"), Vec::new());
    }
    let f = Filaments::new(&all, base.current, SEGS, design::WIRE_RADIUS);
    let tr = tracer(base);
    let Some(ax) = find_axis(&tr, &f, base) else {
        return (done(s, -40.0 - penalty, "no axis"), Vec::new());
    };
    s.axis = Some(ax.rz);
    s.iota_axis = ax.iota(base.n_fp);
    let path = tr.axis_path(&f, ax.rz);

    // Coil points for the clearance test.
    let pts: Vec<[f64; 3]> = all.iter().flat_map(|c| c.sample(SEGS)).collect();
    // Lost: near a coil, or out past the coils (added after run 1: a line escaping outward
    // never nears a coil, and its iota halves could agree at ~0).
    let (centre, out_r) = (tr.centre, tr.lost_radius);
    let lost = |phi: f64, y: Rz| {
        if (y.0 - centre.0).hypot(y.1 - centre.1) > out_r {
            return true;
        }
        let (sp, cp) = phi.sin_cos();
        let p = [y.0 * cp, y.0 * sp, y.1];
        let c2 = CLEARANCE * CLEARANCE;
        pts.iter()
            .any(|q| (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2) < c2)
    };
    let one = |k: usize| {
            let start = k as f64 * LINE_STEP;
            let l = tr.trace_with(&f, (ax.rz.0 + start, ax.rz.1), &path, PERIODS, &lost);
            let halves = l.iota_halves(base.n_fp);
            let good = l.lost_after.is_none() && (halves.0 - halves.1).abs() < IOTA_TOL;
            LineResult {
                start,
                lost_after: l.lost_after,
                iota: l.iota(base.n_fp).0,
                halves,
                good,
            }
    };
    // All lines at once: twelve keep the twelve threads busy (two batches of six, skipping the
    // outer half after a bad inner line, took 14 s against 9.5 s).
    let lines: Vec<LineResult> = (1..=LINES).into_par_iter().map(one).collect();
    let run = lines.iter().take_while(|l| l.good).count();
    s.r_good = run as f64 * LINE_STEP * 1e3;
    s.iota_min = lines[..run]
        .iter()
        .map(|l| l.iota.abs())
        .fold(s.iota_axis, f64::min);
    let value = s.r_good + IOTA_WEIGHT * s.iota_min.min(IOTA_CAP) - penalty;
    (done(s, value, "traced"), lines)
}

/// `opt-probe`: the baseline through the search's score, line by line, to calibrate it.
pub fn probe() {
    let base = design::baseline();
    // The symmetric construction must reproduce the baseline coil set. Rotated copies of an
    // ellipse can start their parameter half a turn on, so compare what matters: the field.
    let built = Filaments::new(&Design::zero().coils(&base), base.current, 128, 0.0);
    let want = Filaments::new(&base.coils(), base.current, 128, 0.0);
    let mut worst = 0.0_f64;
    for (r, phi, z) in [(0.25, 0.1, 0.0), (0.28, 0.9, 0.03), (0.22, 2.0, -0.02), (0.3, 4.0, 0.05), (0.5, 5.5, 0.2)] {
        let (a, b) = (built.b_cyl(r, phi, z), want.b_cyl(r, phi, z));
        let n = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        let e = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt() / n;
        worst = worst.max(e);
    }
    println!("symmetric construction vs baseline coils: field differs by {worst:.1e} (relative, 5 points)");
    let (s, lines) = evaluate_lines(&Design::zero(), &base);
    println!(
        "baseline: value {:.2}, r_good {:.0} mm, iota axis {:.4}, iota min {:.4}, rc_min {:.1} mm, gap {:.1} mm, penalty {:.1}, {} in {:.1} s",
        s.value, s.r_good, s.iota_axis, s.iota_min, s.rc_min, s.gap, s.penalty, s.note, s.secs
    );
    for l in &lines {
        println!(
            "  {:>4.0} mm  iota {:.5}  halves {:.5} / {:.5} (diff {:.1e})  {}{}",
            l.start * 1e3,
            l.iota,
            l.halves.0,
            l.halves.1,
            (l.halves.0 - l.halves.1).abs(),
            if l.good { "good" } else { "bad" },
            l.lost_after.map_or(String::new(), |p| format!(", lost after {p} periods"))
        );
    }
}
