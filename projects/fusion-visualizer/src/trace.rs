//! `trace`: field lines of the shaped equilibrium in the field of 16 TF coils.
//!
//! The toroidal field comes from the TF coils (cfsem oriented circular filaments), the poloidal
//! field from the equilibrium's Psi: B_R = -dPsi/dZ / (2 pi R), B_Z = dPsi/dR / (2 pi R). Lines are
//! integrated in the toroidal angle (RK4, dR/dphi = R B_R / B_phi, dZ/dphi = R B_Z / B_phi) and
//! cut at phi = 0 for Poincare sections. In an axisymmetric field Psi is constant along a line,
//! so the drift of psiN measures how closely lines stay on flux surfaces; q from the traced
//! turns is compared with q = F0 * contour integral of dl / (R |grad Psi|) on the surface.

use std::f64::consts::PI;
use std::fmt::Write as _;
use std::time::Instant;

use cfsem::MU_0;
use rayon::prelude::*;

use crate::contour::{contours, inside, perimeter};
use crate::equilibrium::{Equilibrium, Options, interp, solve, summarize};
use crate::gs::{Grid, Solver};
use crate::shape;
use crate::tf::{TfSet, coil_circular_field, coil_polygon_field};

const B0: f64 = 1.0; // T at R0

fn tf_set() -> TfSet {
    TfSet::for_field(16, 1.05, 0.75, B0, shape::R0)
}

/// Poloidal field from Psi at (R, Z).
fn bpol(grid: &Grid, psi: &[f64], r: f64, z: f64) -> (f64, f64) {
    let (hr, hz) = (0.25 * grid.dr(), 0.25 * grid.dz());
    let dr = (interp(grid, psi, r + hr, z) - interp(grid, psi, r - hr, z)) / (2.0 * hr);
    let dz = (interp(grid, psi, r, z + hz) - interp(grid, psi, r, z - hz)) / (2.0 * hz);
    let k = 2.0 * PI * r;
    (-dz / k, dr / k)
}

struct Line {
    poincare: Vec<(f64, f64)>,
    q: f64,
    /// psiN averaged over every step (all toroidal angles).
    mean_psin: f64,
    dpsin: f64,
    turns: usize,
    lost: bool,
}

/// Trace from (r, 0, z) for `turns` toroidal turns. `ripple` selects the 16-coil field over the
/// axisymmetric F0 / R.
#[allow(clippy::too_many_arguments)]
fn trace(
    grid: &Grid,
    eq: &Equilibrium,
    tf: &TfSet,
    vessel: &[(f64, f64)],
    start: (f64, f64),
    turns: usize,
    per_turn: usize,
    ripple: bool,
) -> Line {
    let f0 = tf.f0();
    let rhs = |phi: f64, r: f64, z: f64| -> (f64, f64) {
        let (br, bz) = bpol(grid, &eq.psi, r, z);
        let (rr, bphi, zz) = if ripple {
            let (a, b, c) = tf.field_cyl(r, phi, z);
            (br + a, b, bz + c)
        } else {
            (br, f0 / r, bz)
        };
        (r * rr / bphi, r * zz / bphi)
    };
    let psin = |r: f64, z: f64| (eq.axis.2 - interp(grid, &eq.psi, r, z)) / (eq.axis.2 - eq.psi_b);
    let h = 2.0 * PI / per_turn as f64;
    let (mut r, mut z) = start;
    let p0 = psin(r, z);
    let mut dpsin = 0.0_f64;
    let theta = |r: f64, z: f64| (z - eq.axis.1).atan2(r - eq.axis.0);
    let mut th_prev = theta(r, z);
    let mut th_total = 0.0;
    // Toroidal angle at each outboard-midplane crossing (theta through a multiple of 2 pi):
    // q from whole poloidal transits only.
    let mut crossings: Vec<f64> = Vec::new();
    let mut psin_sum = 0.0;
    let mut poincare = vec![start];
    for step in 0..turns * per_turn {
        let phi = step as f64 * h;
        let k1 = rhs(phi, r, z);
        let k2 = rhs(phi + 0.5 * h, r + 0.5 * h * k1.0, z + 0.5 * h * k1.1);
        let k3 = rhs(phi + 0.5 * h, r + 0.5 * h * k2.0, z + 0.5 * h * k2.1);
        let k4 = rhs(phi + h, r + h * k3.0, z + h * k3.1);
        r += h / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
        z += h / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
        if !inside(vessel, (r, z)) {
            return Line {
                poincare,
                q: f64::NAN,
                mean_psin: f64::NAN,
                dpsin,
                turns: step / per_turn,
                lost: true,
            };
        }
        psin_sum += psin(r, z);
        let th = theta(r, z);
        let mut d = th - th_prev;
        if d > PI {
            d -= 2.0 * PI;
        } else if d < -PI {
            d += 2.0 * PI;
        }
        let before = th_total;
        th_total += d;
        th_prev = th;
        if (before / (2.0 * PI)).floor() != (th_total / (2.0 * PI)).floor() {
            // Interpolate the crossing within the step.
            let target = (before.max(th_total) / (2.0 * PI)).floor() * 2.0 * PI;
            let t = (target - before) / (th_total - before);
            crossings.push(phi + t * h);
        }
        if (step + 1) % per_turn == 0 {
            poincare.push((r, z));
            dpsin = dpsin.max((psin(r, z) - p0).abs());
        }
    }
    let q = if crossings.len() >= 2 {
        (crossings[crossings.len() - 1] - crossings[0]) / (2.0 * PI * (crossings.len() - 1) as f64)
    } else {
        f64::NAN
    };
    let _ = th_total;
    let mean_psin = psin_sum / (turns * per_turn) as f64;
    Line {
        poincare,
        q,
        mean_psin,
        dpsin,
        turns,
        lost: false,
    }
}

/// q on the flux surface psiN = s from the equilibrium: F0 * contour integral of dl / (R |grad Psi|).
fn q_surface(grid: &Grid, eq: &Equilibrium, f0: f64, s: f64) -> f64 {
    let level = eq.axis.2 - s * (eq.axis.2 - eq.psi_b);
    let Some(c) = contours(grid, &eq.psi, level)
        .into_iter()
        .filter(|c| c.len() > 8 && inside(c, (eq.axis.0, eq.axis.1)))
        .min_by(|a, b| perimeter(a).total_cmp(&perimeter(b)))
    else {
        return f64::NAN;
    };
    let mut sum = 0.0;
    for w in c.windows(2) {
        let (r, z) = (0.5 * (w[0].0 + w[1].0), 0.5 * (w[0].1 + w[1].1));
        let dl = (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1);
        let (br, bz) = bpol(grid, &eq.psi, r, z);
        let grad = 2.0 * PI * r * br.hypot(bz);
        sum += dl / (r * grad);
    }
    f0 * sum
}

/// Outboard-midplane R where psiN = s (bisection from the axis outwards).
fn start_at(grid: &Grid, eq: &Equilibrium, s: f64) -> f64 {
    let psin = |r: f64| (eq.axis.2 - interp(grid, &eq.psi, r, eq.axis.1)) / (eq.axis.2 - eq.psi_b);
    let (mut lo, mut hi) = (eq.axis.0, eq.axis.0 + 0.6);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if psin(mid) < s { lo = mid } else { hi = mid }
    }
    0.5 * (lo + hi)
}

pub fn run() {
    let tf = tf_set();
    println!(
        "TF: {} circular coils, centre R {} m, radius {} m, {:.1} kA each (B0 {B0} T at R0 {} m)",
        tf.n,
        tf.rc,
        tf.radius,
        tf.current / 1e3,
        shape::R0
    );

    println!(
        "1. One TF coil: cfsem's oriented circular filament vs an m-segment polygon of linear filaments"
    );
    let pts = [(1.0, 0.1, 0.2), (0.6, -0.3, 0.5), (1.6, 0.4, -0.1)];
    for m in [64, 256, 1024] {
        let worst = pts
            .iter()
            .map(|&p| {
                let a = coil_circular_field(&tf, p);
                let b = coil_polygon_field(&tf, m, p);
                ((a.0 - b.0).hypot(a.1 - b.1).hypot(a.2 - b.2)) / a.0.hypot(a.1).hypot(a.2)
            })
            .fold(0.0, f64::max);
        println!("   m = {m:>4}: worst relative difference {worst:.2e}");
    }

    println!(
        "2. Ampere: the toroidal average of B_phi is mu0 N I / (2 pi R) exactly; ripple between coils"
    );
    for r in [0.6, 1.0, 1.4] {
        let n_phi = 512;
        let avg: f64 = (0..n_phi)
            .map(|k| tf.field_cyl(r, 2.0 * PI * k as f64 / n_phi as f64, 0.1).1)
            .sum::<f64>()
            / n_phi as f64;
        let want = MU_0 * tf.n as f64 * tf.current / (2.0 * PI * r);
        let on = tf.field_cyl(r, 0.0, 0.0).1;
        let mid = tf.field_cyl(r, PI / tf.n as f64, 0.0).1;
        println!(
            "   R {r}: <B_phi> {avg:.6} T vs {want:.6} T (rel {:.1e}); ripple (on - between)/(on + between) {:.2e}",
            (avg / want - 1.0).abs(),
            (on - mid) / (on + mid)
        );
    }
    // div B by central differences at an off-plane point.
    let (x, y, z, h) = (0.9, 0.3, 0.25, 1e-5);
    let b = |p| tf.field(p);
    let div = (b((x + h, y, z)).0 - b((x - h, y, z)).0 + b((x, y + h, z)).1 - b((x, y - h, z)).1
        + b((x, y, z + h)).2
        - b((x, y, z - h)).2)
        / (2.0 * h);
    let bm = {
        let v = b((x, y, z));
        v.0.hypot(v.1).hypot(v.2)
    };
    println!(
        "   div B at (0.9, 0.3, 0.25): {:.1e} of |B| / 1 m",
        div.abs() / bm
    );

    println!("3. Equilibrium: '{}'", shape::CASES[4].0);
    let grid = shape::grid();
    let solver = Solver::new(&grid, true);
    let m = shape::machine(shape::CASES[4].1, shape::CASES[4].2);
    let prof = shape::profile();
    let eq = solve(&grid, &solver, &m, &prof, &Options::default(), None);
    let s = summarize(&grid, &eq, &prof);
    println!(
        "   axis ({:.4}, {:+.4}), kappa {:.3}, delta {:+.3}, {} iterations",
        eq.axis.0, eq.axis.1, s.kappa, s.delta_upper, eq.iterations
    );

    let turns = 100;
    let per_turn = tf.n * 24;
    println!(
        "4. Field lines from the outboard midplane, {turns} toroidal turns each, RK4 at {per_turn} steps a turn"
    );
    let starts = [0.1, 0.3, 0.5, 0.7, 0.85, 0.95];
    let t0 = Instant::now();
    let jobs: Vec<(f64, bool)> = starts
        .iter()
        .flat_map(|&s| [(s, false), (s, true)])
        .collect();
    let lines: Vec<Line> = jobs
        .par_iter()
        .map(|&(sv, ripple)| {
            let r = start_at(&grid, &eq, sv);
            trace(
                &grid,
                &eq,
                &tf,
                &m.region,
                (r, eq.axis.1),
                turns,
                per_turn,
                ripple,
            )
        })
        .collect();
    println!(
        "   traced {} lines in {:.1} s",
        lines.len(),
        t0.elapsed().as_secs_f64()
    );
    println!(
        "   psiN    q (surface)  q (axisym. TF)  q (16 coils)  max |dpsiN| axisym / 16 coils   16-coil mean psiN, q (surface) there"
    );
    let f0 = tf.f0();
    for (k, &sv) in starts.iter().enumerate() {
        let (a, b) = (&lines[2 * k], &lines[2 * k + 1]);
        println!(
            "   {sv:.2}    {:.4}       {:.4}          {:.4}        {:.1e} / {:.1e}          {:.4}, {:.4}{}",
            q_surface(&grid, &eq, f0, sv),
            a.q,
            b.q,
            a.dpsin,
            b.dpsin,
            b.mean_psin,
            q_surface(&grid, &eq, f0, b.mean_psin),
            if a.lost || b.lost {
                format!("  (lost after {} / {} turns)", a.turns, b.turns)
            } else {
                String::new()
            }
        );
    }

    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out");
    std::fs::create_dir_all(&out).expect("create out/");
    let path = out.join("poincare.svg");
    poincare_svg(&s.lcfs, &m.region, &lines, &path);
    println!(
        "   Poincare plot (16-coil lines, phi = 0): {}",
        path.display()
    );
}

fn poincare_svg(
    lcfs: &[(f64, f64)],
    vessel: &[(f64, f64)],
    lines: &[Line],
    path: &std::path::Path,
) {
    let (r0, r1, z0, z1) = (0.5, 1.5, -0.75, 0.75);
    let px = 600.0;
    let (w, h) = ((r1 - r0) * px, (z1 - z0) * px + 30.0);
    let x = |r: f64| (r - r0) * px;
    let y = |z: f64| (z1 - z) * px + 30.0;
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0}" height="{h:.0}" font-family="sans-serif" font-size="13"><rect width="100%" height="100%" fill="#fff"/><text x="8" y="18" font-weight="bold">Poincare section at phi = 0: 16 TF coils + equilibrium poloidal field</text>"##
    );
    let poly = |pts: &[(f64, f64)], style: &str, svg: &mut String| {
        let d: Vec<String> = pts
            .iter()
            .map(|p| format!("{:.1},{:.1}", x(p.0), y(p.1)))
            .collect();
        let _ = writeln!(svg, r#"<polyline points="{}" {style}/>"#, d.join(" "));
    };
    let mut v = vessel.to_vec();
    v.push(v[0]);
    poly(
        &v,
        r##"fill="none" stroke="#555" stroke-width="2""##,
        &mut svg,
    );
    poly(
        lcfs,
        r##"fill="none" stroke="#c0392b" stroke-width="1.5" stroke-dasharray="5 3""##,
        &mut svg,
    );
    let colors = [
        "#1f77b4", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#ff7f0e",
    ];
    for (k, l) in lines.iter().enumerate().filter(|(k, _)| k % 2 == 1) {
        let c = colors[(k / 2) % colors.len()];
        for p in &l.poincare {
            let _ = writeln!(
                svg,
                r#"<circle cx="{:.1}" cy="{:.1}" r="1.6" fill="{c}"/>"#,
                x(p.0),
                y(p.1)
            );
        }
    }
    let _ = writeln!(svg, "</svg>");
    std::fs::write(path, svg).expect("write SVG");
}
