//! `stray`: magnet sizing for the demonstrator cart, its stray field against the 0.5 mT
//! (5 gauss) pacemaker line, and the rail flux, for all-same and alternating magnet polarity.
//!
//! Magnets are equivalent current loops through cfsem's oriented circular filaments, with their
//! images in the iron rail (z = 0). Soft iron attracts either pole, so polarity is free to choose:
//! all-same adds the dipoles (far field ~ 1/r^3), alternating cancels them.

use std::f64::consts::PI;
use std::fmt::Write as _;

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_density_circular_filament_cartesian_par;

use crate::em::{force_lorentz, images, magnet_loops};

const BR: f64 = 1.3; // T, N42
const LIMIT: f64 = 0.5e-3; // T, the 5 gauss line
const CART_KG: f64 = 0.5; // cart 300 g + payload 200 g
const CORNERS: [(f64, f64); 4] = [
    (-0.05, -0.025),
    (0.05, -0.025),
    (0.05, 0.025),
    (-0.05, 0.025),
];
const LOOPS: usize = 24;

/// A source loop: (radius, current, centre, normal z).
type Src = (f64, f64, (f64, f64, f64), f64);

#[derive(Clone, Copy)]
struct Disc {
    d: f64,
    h: f64,
}

/// Loops of all magnets and their images: (radius, I, centre, normal z).
fn sources(m: Disc, gap: f64, alternating: bool) -> Vec<Src> {
    let mut v = Vec::new();
    for (k, &(x, y)) in CORNERS.iter().enumerate() {
        let sign = if alternating && (k % 2 == 1) {
            -1.0
        } else {
            1.0
        };
        let (loops, i) = magnet_loops(0.5 * m.d, m.h, BR, LOOPS, gap);
        for &(r, z) in &loops {
            v.push((r, sign * i, (x, y, z), 1.0));
            v.push((r, sign * i, (x, y, -z), 1.0)); // image: same current sense
        }
    }
    v
}

fn field(src: &[Src], pts: &[(f64, f64, f64)]) -> Vec<f64> {
    let n = src.len();
    let r: Vec<f64> = src.iter().map(|s| s.0).collect();
    let i: Vec<f64> = src.iter().map(|s| s.1).collect();
    let (lx, ly, lz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
        src.iter().map(|s| s.2.0).collect(),
        src.iter().map(|s| s.2.1).collect(),
        src.iter().map(|s| s.2.2).collect(),
    );
    let nz: Vec<f64> = src.iter().map(|s| s.3).collect();
    let zero = vec![0.0; n];
    let (px, py, pz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
        pts.iter().map(|p| p.0).collect(),
        pts.iter().map(|p| p.1).collect(),
        pts.iter().map(|p| p.2).collect(),
    );
    let (mut bx, mut by, mut bz) = (
        vec![0.0; pts.len()],
        vec![0.0; pts.len()],
        vec![0.0; pts.len()],
    );
    flux_density_circular_filament_cartesian_par(
        (&r, &i),
        (&lx, &ly, &lz),
        (&zero, &zero, &nz),
        &zero,
        (&px, &py, &pz),
        (&mut bx, &mut by, &mut bz),
    )
    .expect("lengths");
    bx.iter()
        .zip(&by)
        .zip(&bz)
        .map(|((a, b), c)| a.hypot(*b).hypot(*c))
        .collect()
}

/// Vertical force on all magnets from all images (the rail's pull on the whole cart), by the
/// Lorentz force on each magnet loop discretized into `seg` points.
fn cart_lift(m: Disc, gap: f64, alternating: bool) -> f64 {
    let src = sources(m, gap, alternating);
    let imgs: Vec<_> = src.iter().filter(|s| s.2.2 > 0.0).copied().collect();
    let n = imgs.len();
    let r: Vec<f64> = imgs.iter().map(|s| s.0).collect();
    let i: Vec<f64> = imgs.iter().map(|s| s.1).collect();
    let (lx, ly, lz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
        imgs.iter().map(|s| s.2.0).collect(),
        imgs.iter().map(|s| s.2.1).collect(),
        imgs.iter().map(|s| s.2.2).collect(),
    );
    let zero = vec![0.0; n];
    let one = vec![1.0; n];
    let seg = 96;
    let mut f = 0.0;
    for s in src.iter().filter(|s| s.2.2 < 0.0) {
        let (rad, cur, (cx, cy, cz), _) = *s;
        let pts: Vec<(f64, f64, f64)> = (0..seg)
            .map(|k| {
                let t = 2.0 * PI * (k as f64 + 0.5) / seg as f64;
                (cx + rad * t.cos(), cy + rad * t.sin(), cz)
            })
            .collect();
        let (px, py, pz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
            pts.iter().map(|p| p.0).collect(),
            pts.iter().map(|p| p.1).collect(),
            pts.iter().map(|p| p.2).collect(),
        );
        let (mut bx, mut by, mut bz) = (vec![0.0; seg], vec![0.0; seg], vec![0.0; seg]);
        flux_density_circular_filament_cartesian_par(
            (&r, &i),
            (&lx, &ly, &lz),
            (&zero, &zero, &one),
            &zero,
            (&px, &py, &pz),
            (&mut bx, &mut by, &mut bz),
        )
        .expect("lengths");
        // dl = (-sin t, cos t, 0) rad dt; F_z = I (dl x B)_z = I (dl_x B_y - dl_y B_x).
        for k in 0..seg {
            let t = 2.0 * PI * (k as f64 + 0.5) / seg as f64;
            let dl = 2.0 * PI * rad / seg as f64;
            f += cur * dl * (-t.sin() * by[k] - t.cos() * bx[k]);
        }
        let _ = bz;
    }
    f
}

/// Gap (m) at which one magnet alone carries `need` (N), by bisection; None if it cannot at 1 mm.
fn gap_for(m: Disc, need: f64) -> Option<f64> {
    let lift = |g: f64| {
        let (l, i) = magnet_loops(0.5 * m.d, m.h, BR, LOOPS, g);
        force_lorentz(&l, i, &images(&l), i)
    };
    if lift(0.001) < need {
        return None;
    }
    let (mut lo, mut hi) = (0.001, 0.05);
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        if lift(mid) > need { lo = mid } else { hi = mid }
    }
    Some(0.5 * (lo + hi))
}

/// Distance from `start` along unit `dir` at which |B| falls to the limit (m).
fn reach(src: &[Src], start: (f64, f64, f64), dir: (f64, f64, f64)) -> f64 {
    let at = |s: f64| {
        field(
            src,
            &[(
                start.0 + s * dir.0,
                start.1 + s * dir.1,
                start.2 + s * dir.2,
            )],
        )[0]
    };
    let mut hi = 0.01;
    while at(hi) > LIMIT && hi < 5.0 {
        hi *= 1.5;
    }
    let mut lo = 0.0;
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        if at(mid) > LIMIT { lo = mid } else { hi = mid }
    }
    0.5 * (lo + hi)
}

pub fn run() {
    let need = CART_KG * 9.81 / 4.0;
    println!(
        "Demonstrator cart: {CART_KG} kg on 4 N42 discs (Br {BR} T) at the corners of 100 x 50 mm: \
         {need:.3} N a corner"
    );
    println!("1. Gap at which one disc carries its corner (stock disc sizes, diameter x height)");
    let discs = [
        Disc { d: 0.010, h: 0.003 },
        Disc { d: 0.012, h: 0.004 },
        Disc { d: 0.015, h: 0.005 },
        Disc { d: 0.020, h: 0.005 },
    ];
    let mut chosen = None;
    for m in discs {
        let g = gap_for(m, need);
        let moment = BR / MU_0 * PI * (0.5 * m.d).powi(2) * m.h;
        match g {
            Some(g) => println!(
                "   {:>2.0} x {:.0} mm: gap {:.2} mm, dipole moment {moment:.3} A m^2",
                m.d * 1e3,
                m.h * 1e3,
                g * 1e3
            ),
            None => println!(
                "   {:>2.0} x {:.0} mm: cannot carry it even at 1 mm",
                m.d * 1e3,
                m.h * 1e3
            ),
        }
        // The demonstrator: the smallest disc that floats the corner at 5 mm or more.
        if chosen.is_none() && g.is_some_and(|g| g >= 0.005) {
            chosen = Some((m, g.expect("checked")));
        }
    }
    let (m, gap) = chosen.expect("a disc that floats at 5 mm");
    println!(
        "   chosen: {:.0} x {:.0} mm at {:.2} mm",
        m.d * 1e3,
        m.h * 1e3,
        gap * 1e3
    );

    println!("2. Whole cart (each magnet against every image) and the 0.5 mT reach");
    println!(
        "   layout       lift (4 corners)  vs 4 x single  0.5 mT reach: down from centre / down from a magnet / sideways from a magnet's edge / 45 deg"
    );
    let mut maps = Vec::new();
    for alternating in [false, true] {
        let src = sources(m, gap, alternating);
        let lift = cart_lift(m, gap, alternating);
        let zc = -gap - m.h; // magnet bottom faces
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let r_down = reach(&src, (0.0, 0.0, zc), (0.0, 0.0, -1.0));
        let r_mag = reach(&src, (0.05, 0.025, zc), (0.0, 0.0, -1.0));
        let r_side = reach(
            &src,
            (0.05, 0.025 + 0.5 * m.d, -gap - 0.5 * m.h),
            (0.0, 1.0, 0.0),
        );
        let r_diag = reach(&src, (0.05, 0.025, zc), (s, 0.0, -s));
        println!(
            "   {:<11}  {lift:.3} N          {:+.2}%        {:.0} / {:.0} / {:.0} / {:.0} mm",
            if alternating {
                "alternating"
            } else {
                "all same"
            },
            (lift / (4.0 * need) - 1.0) * 100.0,
            r_down * 1e3,
            r_mag * 1e3,
            r_side * 1e3,
            r_diag * 1e3
        );
        maps.push((alternating, src));
    }

    println!("3. Rail flux: total flux entering the rail from one magnet (alternating layout)");
    {
        let src = &maps[1].1;
        // B_z just below the plane over the quadrant around the (0.05, 0.025) magnet's sign region.
        let n = 160;
        let (x0, x1, y0, y1) = (-0.15, 0.15, -0.1, 0.1);
        let (dx, dy) = ((x1 - x0) / n as f64, (y1 - y0) / n as f64);
        let pts: Vec<(f64, f64, f64)> = (0..n)
            .flat_map(|a| {
                (0..n).map(move |b| {
                    (
                        x0 + (a as f64 + 0.5) * dx,
                        y0 + (b as f64 + 0.5) * dy,
                        -1e-6,
                    )
                })
            })
            .collect();
        let bz = field_z(src, &pts);
        let up: f64 = bz.iter().filter(|b| **b > 0.0).sum::<f64>() * dx * dy;
        let per = up / 2.0; // two magnets share each sign
        let width = 0.06;
        let t = per / (1.2 * width);
        println!(
            "   {:.3e} Wb of each sign through the rail face, {:.3e} Wb a magnet; in a {:.0} mm wide \
             rail kept below 1.2 T that needs >= {:.2} mm of steel",
            up,
            per,
            width * 1e3,
            t * 1e3
        );
    }

    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out");
    std::fs::create_dir_all(&out).expect("create out/");
    for (alternating, src) in &maps {
        let name = if *alternating {
            "stray_alternating.svg"
        } else {
            "stray_same.svg"
        };
        map_svg(src, m, gap, *alternating, &out.join(name));
    }
    println!(
        "4. Maps (x-z plane through a row of magnets): {}",
        out.join("stray_*.svg").display()
    );
}

fn field_z(src: &[Src], pts: &[(f64, f64, f64)]) -> Vec<f64> {
    let n = src.len();
    let r: Vec<f64> = src.iter().map(|s| s.0).collect();
    let i: Vec<f64> = src.iter().map(|s| s.1).collect();
    let (lx, ly, lz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
        src.iter().map(|s| s.2.0).collect(),
        src.iter().map(|s| s.2.1).collect(),
        src.iter().map(|s| s.2.2).collect(),
    );
    let zero = vec![0.0; n];
    let one = vec![1.0; n];
    let (px, py, pz): (Vec<f64>, Vec<f64>, Vec<f64>) = (
        pts.iter().map(|p| p.0).collect(),
        pts.iter().map(|p| p.1).collect(),
        pts.iter().map(|p| p.2).collect(),
    );
    let (mut bx, mut by, mut bz) = (
        vec![0.0; pts.len()],
        vec![0.0; pts.len()],
        vec![0.0; pts.len()],
    );
    flux_density_circular_filament_cartesian_par(
        (&r, &i),
        (&lx, &ly, &lz),
        (&zero, &zero, &one),
        &zero,
        (&px, &py, &pz),
        (&mut bx, &mut by, &mut bz),
    )
    .expect("lengths");
    bz
}

/// |B| on the x-z plane through the y = 0.025 row, with 0.5 / 1 / 5 mT contours.
fn map_svg(src: &[Src], m: Disc, gap: f64, alternating: bool, path: &std::path::Path) {
    let (nx, nz) = (241, 151);
    let (x0, x1, z0, z1) = (-0.15, 0.15, -0.12, 0.0);
    let pts: Vec<(f64, f64, f64)> = (0..nx)
        .flat_map(|a| {
            (0..nz).map(move |b| {
                (
                    x0 + (x1 - x0) * a as f64 / (nx - 1) as f64,
                    0.025,
                    z0 + (z1 - z0) * b as f64 / (nz - 1) as f64 - 1e-4,
                )
            })
        })
        .collect();
    let b = field(src, &pts);
    let px = 3000.0;
    let (w, h) = ((x1 - x0) * px, (z1 - z0) * px + 40.0);
    let sx = |x: f64| (x - x0) * px;
    let sz = |z: f64| (z1 - z) * px + 40.0;
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0}" height="{h:.0}" font-family="sans-serif" font-size="13"><rect width="100%" height="100%" fill="#fff"/><text x="8" y="18" font-weight="bold">|B| below the rail, {} polarity: {:.0} x {:.0} mm N42 discs at {:.1} mm; contours 5 / 1 / 0.5 mT</text>"##,
        if alternating {
            "alternating"
        } else {
            "all-same"
        },
        m.d * 1e3,
        m.h * 1e3,
        gap * 1e3
    );
    // Marching squares on each level.
    let at = |a: usize, c: usize| b[a * nz + c];
    let pos = |a: usize, c: usize| {
        (
            x0 + (x1 - x0) * a as f64 / (nx - 1) as f64,
            z0 + (z1 - z0) * c as f64 / (nz - 1) as f64,
        )
    };
    for (level, color, width) in [
        (5e-3, "#c0392b", 1.5),
        (1e-3, "#e67e22", 1.5),
        (0.5e-3, "#2e86c1", 2.5),
    ] {
        for a in 0..nx - 1 {
            for c in 0..nz - 1 {
                let corners = [(a, c), (a + 1, c), (a + 1, c + 1), (a, c + 1)];
                let v: Vec<f64> = corners.iter().map(|&(p, q)| at(p, q) - level).collect();
                let mut cross = Vec::new();
                for e in 0..4 {
                    let (i, j) = (e, (e + 1) % 4);
                    if (v[i] < 0.0) != (v[j] < 0.0) {
                        let t = v[i] / (v[i] - v[j]);
                        let (p, q) = (
                            pos(corners[i].0, corners[i].1),
                            pos(corners[j].0, corners[j].1),
                        );
                        cross.push((p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1)));
                    }
                }
                for pair in cross.chunks(2).filter(|p| p.len() == 2) {
                    let _ = writeln!(
                        svg,
                        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{color}" stroke-width="{width}"/>"#,
                        sx(pair[0].0),
                        sz(pair[0].1),
                        sx(pair[1].0),
                        sz(pair[1].1)
                    );
                }
            }
        }
    }
    // Rail, magnets (shaded by polarity).
    let _ = writeln!(
        svg,
        r##"<rect x="0" y="{:.1}" width="{w:.0}" height="12" fill="#666"/><text x="8" y="{:.1}" fill="#fff">steel rail</text>"##,
        sz(0.0) - 12.0,
        sz(0.0) - 2.0
    );
    for (k, &(x, y)) in CORNERS.iter().enumerate() {
        if (y - 0.025).abs() > 1e-9 {
            continue;
        }
        let up = !(alternating && k % 2 == 1);
        let _ = writeln!(
            svg,
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}"/>"##,
            sx(x - 0.5 * m.d),
            sz(-gap),
            m.d * px,
            m.h * px,
            if up { "#c0392b" } else { "#2e86c1" }
        );
    }
    let _ = writeln!(svg, "</svg>");
    std::fs::write(path, svg).expect("write SVG");
}
