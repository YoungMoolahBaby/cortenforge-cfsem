//! Plain SVG plots of an equilibrium in the (R, Z) plane.

use std::fmt::Write as _;

use crate::contour::contours;
use crate::equilibrium::{Boundary, Equilibrium, Machine, Summary};
use crate::gs::Grid;

const PX_PER_M: f64 = 400.0;

pub fn equilibrium_svg(
    grid: &Grid,
    eq: &Equilibrium,
    s: &Summary,
    m: &Machine,
    title: &str,
    path: &std::path::Path,
) {
    // View: the grid plus the coils.
    let mut r0 = grid.rs[0];
    let mut r1 = *grid.rs.last().expect("grid");
    let mut z0 = grid.zs[0];
    let mut z1 = *grid.zs.last().expect("grid");
    for c in &m.coils {
        r0 = r0.min(c.r - c.w);
        r1 = r1.max(c.r + c.w);
        z0 = z0.min(c.z - c.h);
        z1 = z1.max(c.z + c.h);
    }
    let pad = 0.05;
    let (r0, r1, z0, z1) = (r0.max(0.0) - pad, r1 + pad, z0 - pad, z1 + pad);
    let (w, h) = ((r1 - r0) * PX_PER_M, (z1 - z0) * PX_PER_M + 60.0);
    let x = |r: f64| (r - r0) * PX_PER_M;
    let y = |z: f64| (z1 - z) * PX_PER_M + 40.0;
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0}" height="{h:.0}" viewBox="0 0 {w:.0} {h:.0}" font-family="sans-serif" font-size="13">"#
    );
    let _ = writeln!(
        svg,
        r##"<rect width="100%" height="100%" fill="#ffffff"/>"##
    );
    let _ = writeln!(
        svg,
        r##"<text x="8" y="18" font-weight="bold">{title}</text>"##
    );
    let _ = writeln!(
        svg,
        r##"<text x="8" y="34" fill="#444">R_geo {:.3} m, a {:.3} m, kappa {:.3}, delta {:+.3}/{:+.3}, beta_p {:.3}, li {:.3}, Ip {:.0} kA, {}</text>"##,
        s.r_geo,
        s.a,
        s.kappa,
        s.delta_upper,
        s.delta_lower,
        s.beta_p,
        s.li,
        s.ip / 1e3,
        match eq.boundary {
            Boundary::Limiter(..) => "limited",
            Boundary::XPoint(..) => "diverted",
        }
    );
    // Grid outline.
    let _ = writeln!(
        svg,
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#bbb" stroke-dasharray="4 3"/>"##,
        x(grid.rs[0]),
        y(*grid.zs.last().expect("grid")),
        x(*grid.rs.last().expect("grid")) - x(grid.rs[0]),
        y(grid.zs[0]) - y(*grid.zs.last().expect("grid"))
    );
    let poly = |pts: &[(f64, f64)], style: &str, svg: &mut String| {
        let d: Vec<String> = pts
            .iter()
            .map(|p| format!("{:.1},{:.1}", x(p.0), y(p.1)))
            .collect();
        let _ = writeln!(svg, r#"<polyline points="{}" {style}/>"#, d.join(" "));
    };
    // Flux contours: inside the plasma evenly in psiN, outside a few more.
    let span = eq.axis.2 - eq.psi_b;
    for k in 1..20 {
        let psin = k as f64 * 0.1;
        let level = eq.axis.2 - psin * span;
        let style = if psin < 1.0 {
            r##"fill="none" stroke="#3b6fb6" stroke-width="1""##
        } else {
            r##"fill="none" stroke="#9aa" stroke-width="0.8" stroke-dasharray="3 3""##
        };
        for c in contours(grid, &eq.psi, level) {
            // Plasma levels only inside the last closed surface (coils reach the same flux).
            if psin < 1.0 && !crate::contour::inside(&s.lcfs, c[0]) {
                continue;
            }
            poly(&c, style, &mut svg);
        }
    }
    // Vessel, limiter points, last closed surface.
    let mut vessel = m.region.clone();
    vessel.push(vessel[0]);
    poly(
        &vessel,
        r##"fill="none" stroke="#555" stroke-width="2""##,
        &mut svg,
    );
    poly(
        &s.lcfs,
        r##"fill="none" stroke="#c0392b" stroke-width="2.5""##,
        &mut svg,
    );
    // Coils, shaded by current sign.
    for c in &m.coils {
        let fill = if c.current >= 0.0 {
            "#e67e22"
        } else {
            "#2e86c1"
        };
        let _ = writeln!(
            svg,
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{fill}" stroke="#333"/><text x="{:.1}" y="{:.1}" font-size="11">{} {:+.0} kA</text>"##,
            x(c.r - 0.5 * c.w),
            y(c.z + 0.5 * c.h),
            c.w * PX_PER_M,
            c.h * PX_PER_M,
            x(c.r + 0.5 * c.w) + 3.0,
            y(c.z) + 4.0,
            c.name,
            c.current / 1e3
        );
    }
    // Axis and X-point.
    let _ = writeln!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" fill="#c0392b" text-anchor="middle">+</text>"##,
        x(eq.axis.0),
        y(eq.axis.1) + 5.0
    );
    if let Boundary::XPoint(r, z) = eq.boundary {
        let _ = writeln!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" fill="#c0392b" text-anchor="middle" font-weight="bold">x</text>"##,
            x(r),
            y(z) + 5.0
        );
    }
    let _ = writeln!(svg, "</svg>");
    std::fs::write(path, svg).expect("write SVG");
}
