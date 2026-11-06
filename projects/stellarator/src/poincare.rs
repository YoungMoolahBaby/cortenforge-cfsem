//! `trace`: the baseline's magnetic axis, field lines from the axis outward, the rotational
//! transform profile, and Poincare sections at phi = 0 and the half-period plane (SVG).

use std::fmt::Write as _;
use std::time::Instant;

use rayon::prelude::*;

use crate::coils::{Filaments, FourierCurve, RotatingEllipse};
use crate::design;
use crate::search;
use crate::trace::{Line, Rz};

/// Field periods traced per line.
const PERIODS: usize = 600;

pub fn run() {
    let c = design::baseline();
    trace_coils(&c, &c.coils(), "baseline", "poincare.svg");
}

/// Axis, lines, iota profile and Poincare sections for any coil set built on `c`'s layout
/// (field periods, centre, current). A line is lost when it comes within
/// `search::CLEARANCE` of a coil centreline, the same rule the coil search uses.
pub fn trace_coils(c: &RotatingEllipse, coils: &[FourierCurve], label: &str, file: &str) {
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out");
    std::fs::create_dir_all(&out).expect("create out");
    let f = Filaments::new(coils, c.current, design::SEGS, design::WIRE_RADIUS);
    let tr = design::tracer(c);
    let pts: Vec<[f64; 3]> = coils.iter().flat_map(|k| k.sample(design::SEGS)).collect();
    let c2 = search::CLEARANCE * search::CLEARANCE;
    // Lost: near a coil, or out past the coils (a line escaping outward never nears one).
    let (centre, out_r) = (tr.centre, tr.lost_radius);
    let lost = |phi: f64, y: Rz| {
        if (y.0 - centre.0).hypot(y.1 - centre.1) > out_r {
            return true;
        }
        let (sp, cp) = phi.sin_cos();
        let p = [y.0 * cp, y.0 * sp, y.1];
        pts.iter()
            .any(|q| (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2) < c2)
    };
    println!("{label}:");
    println!(
        "  layout: n_fp {}, {} coils, rc {} m, starting ellipse {} x {} m, helix shift R {:+} Z {:+} m, {:.2} kA-turns, B0 {} T",
        c.n_fp,
        c.n_coils(),
        c.rc,
        c.a_long,
        c.a_short,
        c.shift,
        c.shift_z,
        c.current / 1e3,
        design::B0
    );

    let t0 = Instant::now();
    let fps = tr.fixed_points(&f, (c.rc, 0.0));
    for p in &fps {
        println!(
            "  fixed point R {:.5} Z {:+.2e}: {} (trace {:+.6}, det {:.9})",
            p.rz.0,
            p.rz.1,
            if p.elliptic() { "elliptic" } else { "hyperbolic" },
            p.trace,
            p.det
        );
    }
    let Some(ax) = fps.iter().find(|p| p.elliptic()).copied() else {
        println!("  no elliptic axis");
        return;
    };
    let path = tr.axis_path(&f, ax.rz);
    let half = path[tr.steps / 2 - 1];
    println!(
        "axis at phi = 0: R {:.5} m; at the half period: R {:.5} Z {:+.5} m; iota on axis {:.5} ({:.1} s)",
        ax.rz.0,
        half.0,
        half.1,
        ax.iota(c.n_fp),
        t0.elapsed().as_secs_f64()
    );

    // Lines from the axis outward along the phi = 0 midplane.
    let t0 = Instant::now();
    let starts: Vec<f64> = (1..=36).map(|k| 0.0025 * k as f64).collect();
    let lines: Vec<Line> = starts
        .par_iter()
        .map(|d| tr.trace_with(&f, (ax.rz.0 + d, 0.0), &path, PERIODS, &lost))
        .collect();
    println!(
        "{} lines x {} periods ({:.0} toroidal turns) in {:.1} s",
        lines.len(),
        PERIODS,
        PERIODS as f64 / c.n_fp as f64,
        t0.elapsed().as_secs_f64()
    );
    println!("  start (mm)   iota (Birkhoff)   iota (mean)   rho min-max (mm)   height/width phi=0, half   fate");
    let mut last_confined = None;
    for (d, l) in starts.iter().zip(&lines) {
        let (iota, mean) = l.iota(c.n_fp);
        let fate = match l.lost_after {
            Some(p) => format!("lost after {p} periods"),
            None => {
                last_confined = Some(*d);
                "confined".into()
            }
        };
        println!(
            "  {:>9.1}   {:>15.5}   {:>11.5}   {:>7.1} - {:<7.1}   {:>8.3} {:>8.3}          {fate}",
            d * 1e3,
            iota,
            mean,
            l.rho_min * 1e3,
            l.rho_max * 1e3,
            aspect(&l.plane0),
            aspect(&l.plane_half)
        );
    }
    if let Some(d) = last_confined {
        println!("outermost confined start: {:.1} mm from the axis", d * 1e3);
    }

    let svg = svg(c, label, &ax.rz, half, &starts, &lines);
    let p = out.join(file);
    std::fs::write(&p, svg).expect("write svg");
    println!("wrote {}", p.display());
}

/// Two Poincare panels (phi = 0, half period) and the iota profile.
fn svg(
    c: &RotatingEllipse,
    label: &str,
    axis: &Rz,
    half: Rz,
    starts: &[f64],
    lines: &[Line],
) -> String {
    let px = 2600.0; // pixels per metre
    let span = 2.0 * (c.a_long + c.shift.abs()) + 0.02;
    let panel = span * px;
    let (w, h) = (3.0 * panel + 80.0, panel + 90.0);
    let mut s = String::new();
    let _ = writeln!(
        s,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0}" height="{h:.0}" viewBox="0 0 {w:.0} {h:.0}" font-family="sans-serif" font-size="13">"#
    );
    let _ = writeln!(s, r##"<rect width="100%" height="100%" fill="#ffffff"/>"##);
    let _ = writeln!(
        s,
        r##"<text x="10" y="20" font-weight="bold">Toy stellarator ({label}): n_fp {}, {} modular coils, Poincare sections ({} periods a line) and rotational transform</text>"##,
        c.n_fp,
        c.n_coils(),
        PERIODS
    );
    let colour = |k: usize, lost: bool| {
        if lost {
            "#bbbbbb".to_string()
        } else {
            let hue = 240.0 - 240.0 * k as f64 / starts.len().max(1) as f64;
            format!("hsl({hue:.0},75%,42%)")
        }
    };
    for (panel_i, (title, ax)) in [("phi = 0", *axis), ("phi = half period", half)].iter().enumerate() {
        let x0 = 20.0 + panel_i as f64 * (panel + 20.0);
        let y0 = 50.0;
        let (rc, zc) = (c.rc, 0.0);
        let x = |r: f64| x0 + (r - rc + 0.5 * span) * px;
        let y = |z: f64| y0 + (zc - z + 0.5 * span) * px;
        let _ = writeln!(
            s,
            r##"<rect x="{x0:.1}" y="{y0:.1}" width="{panel:.1}" height="{panel:.1}" fill="none" stroke="#999"/>"##
        );
        let _ = writeln!(s, r##"<text x="{:.1}" y="{:.1}">{title}</text>"##, x0 + 6.0, y0 + 16.0);
        for (k, l) in lines.iter().enumerate() {
            let col = colour(k, l.lost_after.is_some());
            let pts = if panel_i == 0 { &l.plane0 } else { &l.plane_half };
            for p in pts {
                let (px_, py_) = (x(p.0), y(p.1));
                if px_ > x0 && px_ < x0 + panel && py_ > y0 && py_ < y0 + panel {
                    let _ = writeln!(s, r#"<circle cx="{px_:.1}" cy="{py_:.1}" r="0.9" fill="{col}"/>"#);
                }
            }
        }
        let _ = writeln!(
            s,
            r##"<path d="M{:.1} {:.1}h-6h12m-6 -6v12" stroke="#000" stroke-width="1.5"/>"##,
            x(ax.0),
            y(ax.1)
        );
        // 2 cm scale bar.
        let _ = writeln!(
            s,
            r##"<path d="M{:.1} {:.1}h{:.1}" stroke="#000" stroke-width="2"/><text x="{:.1}" y="{:.1}">2 cm</text>"##,
            x0 + 10.0,
            y0 + panel - 12.0,
            0.02 * px,
            x0 + 14.0 + 0.02 * px,
            y0 + panel - 8.0
        );
    }
    // Iota profile against the start distance.
    let x0 = 20.0 + 2.0 * (panel + 20.0);
    let y0 = 50.0;
    let dmax = starts.last().copied().unwrap_or(0.05) * 1.05;
    let imax = lines
        .iter()
        .filter(|l| l.lost_after.is_none())
        .map(|l| l.iota(c.n_fp).0.abs())
        .fold(0.0, f64::max)
        * 1.15;
    let x = |d: f64| x0 + 40.0 + d / dmax * (panel - 50.0);
    let y = |i: f64| y0 + panel - 30.0 - i / imax * (panel - 50.0);
    let _ = writeln!(
        s,
        r##"<rect x="{x0:.1}" y="{y0:.1}" width="{panel:.1}" height="{panel:.1}" fill="none" stroke="#999"/>"##
    );
    let _ = writeln!(
        s,
        r##"<text x="{:.1}" y="{:.1}">|iota| vs start distance from the axis (grey: lost)</text>"##,
        x0 + 6.0,
        y0 + 16.0
    );
    let _ = writeln!(
        s,
        r##"<path d="M{:.1} {:.1}V{:.1}H{:.1}" fill="none" stroke="#000"/>"##,
        x(0.0),
        y(imax),
        y(0.0),
        x(dmax)
    );
    for k in 0..=5 {
        let iv = imax * k as f64 / 5.0;
        let _ = writeln!(s, r##"<text x="{:.1}" y="{:.1}" text-anchor="end">{iv:.3}</text>"##, x(0.0) - 4.0, y(iv) + 4.0);
    }
    for k in 0..=5 {
        let dv = dmax * k as f64 / 5.0;
        let _ = writeln!(s, r##"<text x="{:.1}" y="{:.1}" text-anchor="middle">{:.0} mm</text>"##, x(dv), y(0.0) + 16.0, dv * 1e3);
    }
    for (k, (d, l)) in starts.iter().zip(lines).enumerate() {
        let col = colour(k, l.lost_after.is_some());
        let i = l.iota(c.n_fp).0.abs();
        if i.is_finite() {
            let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{:.1}" r="3.5" fill="{col}"/>"#, x(*d), y(i.min(imax)));
        }
    }
    let _ = writeln!(s, "</svg>");
    s
}

/// Height over width of a set of section points (Z extent / R extent).
fn aspect(pts: &[Rz]) -> f64 {
    let (mut r0, mut r1, mut z0, mut z1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for p in pts {
        r0 = r0.min(p.0);
        r1 = r1.max(p.0);
        z0 = z0.min(p.1);
        z1 = z1.max(p.1);
    }
    (z1 - z0) / (r1 - r0)
}
