//! `mesh`: the baseline's coils swept into solids with cf-design's `Solid::pipe_spline`,
//! meshed, checked (watertight, manifold, volume against pi r^2 L) and written as STL (mm).
//!
//! `pipe_spline` takes an open Catmull-Rom curve (end tangents clamped). A closed coil goes in
//! with its control points wrapped: [p_(n-1), p_0, ..., p_(n-1), p_0, p_1]. Every true span then
//! appears with its real neighbours; the two extra end spans repeat true spans with clamped
//! tangents.

use std::f64::consts::PI;
use std::time::Instant;

use cortenforge::cf_design::Solid;
use cortenforge::mesh::io::save_mesh;
use cortenforge::mesh::repair::{remove_unreferenced_vertices, validate_mesh, weld_vertices};
use cortenforge::mesh::types::{IndexedMesh, Point3};
use rayon::prelude::*;

use crate::coils::{FourierCurve, min_distance};
use crate::design;
use crate::render;

const MM: f64 = 1e3;

fn control_points(c: &FourierCurve, n: usize, wrap: bool) -> Vec<Point3<f64>> {
    let p: Vec<Point3<f64>> = c
        .sample(n)
        .into_iter()
        .map(|q| Point3::new(q[0] * MM, q[1] * MM, q[2] * MM))
        .collect();
    if wrap {
        let mut v = vec![p[n - 1]];
        v.extend_from_slice(&p);
        v.push(p[0]);
        v.push(p[1]);
        v
    } else {
        let mut v = p.clone();
        v.push(p[0]);
        v
    }
}

fn merge(parts: &[IndexedMesh]) -> IndexedMesh {
    let mut out = IndexedMesh::new();
    for p in parts {
        let base = out.vertices.len() as u32;
        out.vertices.extend_from_slice(&p.vertices);
        out.faces.extend(p.faces.iter().map(|f| f.map(|i| i + base)));
    }
    out
}

struct Part {
    mesh: IndexedMesh,
    secs: f64,
}

/// Move every vertex onto the solid's zero level set along its SDF gradient (central
/// differences), `iters` Newton steps. The octree mesher's topology is kept.
fn snap(s: &Solid, m: &mut IndexedMesh, iters: usize) {
    let h = 1e-3;
    m.vertices.par_iter_mut().for_each(|v| {
        for _ in 0..iters {
            let d = s.evaluate(v);
            let g = |dx: f64, dy: f64, dz: f64| {
                (s.evaluate(&Point3::new(v.x + dx, v.y + dy, v.z + dz))
                    - s.evaluate(&Point3::new(v.x - dx, v.y - dy, v.z - dz)))
                    / (2.0 * h)
            };
            let (gx, gy, gz) = (g(h, 0.0, 0.0), g(0.0, h, 0.0), g(0.0, 0.0, h));
            let n2 = gx * gx + gy * gy + gz * gz;
            if n2 > 1e-12 {
                *v = Point3::new(v.x - d * gx / n2, v.y - d * gy / n2, v.z - d * gz / n2);
            }
        }
    });
}

/// Sweep a coil with `pipe_spline` (wrapped control points), mesh it with the octree mesher and
/// optionally snap the vertices onto the surface. No weld: the octree mesh is already clean, and
/// welding at 1 um made it non-manifold (see `mesh-probe`). The uniform marching-cubes `mesh` took
/// 305 s for one coil at 1 mm, against 1.7 s here.
fn sweep(c: &FourierCurve, n: usize, r_mm: f64, cell: f64, snapped: bool) -> Part {
    let t0 = Instant::now();
    let s = Solid::pipe_spline(control_points(c, n, true), r_mm);
    let mut mesh = s.mesh_adaptive_par(cell).geometry;
    if snapped {
        snap(&s, &mut mesh, 2);
    }
    Part {
        mesh,
        secs: t0.elapsed().as_secs_f64(),
    }
}

fn report(name: &str, p: &Part, exact: f64, curve: &FourierCurve, r_mm: f64) {
    let r = validate_mesh(&p.mesh);
    let vol = p.mesh.signed_volume();
    let (lo, hi) = deviation(&p.mesh, curve, r_mm);
    println!(
        "  {name:<34} {:>7} faces  watertight {:<5} manifold {:<5}  volume vs pi r^2 L {:+.2e}  surface {:+.3} to {:+.3} mm  {:.1} s",
        p.mesh.faces.len(),
        r.is_watertight,
        r.is_manifold,
        vol / exact - 1.0,
        lo,
        hi,
        p.secs
    );
}

/// Range of (distance from the true centreline - r) over the mesh vertices, in mm.
fn deviation(m: &IndexedMesh, curve: &FourierCurve, r_mm: f64) -> (f64, f64) {
    let line: Vec<[f64; 3]> = curve.sample(8192).into_iter().map(|q| q.map(|x| x * MM)).collect();
    let d: Vec<f64> = m
        .vertices
        .par_iter()
        .map(|v| {
            let mut d2 = f64::INFINITY;
            for q in &line {
                d2 = d2.min((v.x - q[0]).powi(2) + (v.y - q[1]).powi(2) + (v.z - q[2]).powi(2));
            }
            d2.sqrt() - r_mm
        })
        .collect();
    (d.iter().copied().fold(f64::MAX, f64::min), d.iter().copied().fold(f64::MIN, f64::max))
}

pub fn run() {
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out").join("mesh");
    std::fs::create_dir_all(&out).expect("create out/mesh");
    let c = design::baseline();
    let coils = c.coils();
    let r_mm = design::WIRE_RADIUS * MM;
    let exact = |k: usize| PI * r_mm * r_mm * coils[k].length(256) * MM;

    println!("probe: coil 1 ({:.1} mm long), conductor radius {r_mm} mm", coils[1].length(256) * MM);
    for (n, cell, snapped) in [(64, 1.0, false), (128, 1.0, false), (128, 1.0, true), (128, 0.5, true)] {
        let p = sweep(&coils[1], n, r_mm, cell, snapped);
        report(
            &format!("{n} points, cell {cell} mm{}", if snapped { ", snapped" } else { "" }),
            &p,
            exact(1),
            &coils[1],
            r_mm,
        );
    }

    // All coils at the chosen settings, one at a time (the mesher is parallel).
    let (n, cell) = (128, 1.0);
    let t0 = Instant::now();
    let parts: Vec<Part> = coils.iter().map(|cc| sweep(cc, n, r_mm, cell, true)).collect();
    println!(
        "all {} coils ({n} points, cell {cell} mm, snapped) in {:.1} s:",
        coils.len(),
        t0.elapsed().as_secs_f64()
    );
    let (mut worst_v, mut worst_s) = (0.0_f64, 0.0_f64);
    let mut all_ok = true;
    for (k, p) in parts.iter().enumerate() {
        let r = validate_mesh(&p.mesh);
        all_ok &= r.is_watertight && r.is_manifold;
        worst_v = worst_v.max((p.mesh.signed_volume() / exact(k) - 1.0).abs());
        let (lo, hi) = deviation(&p.mesh, &coils[k], r_mm);
        worst_s = worst_s.max(lo.abs()).max(hi.abs());
        if k < c.coils_per_period {
            report(&format!("coil {k}"), p, exact(k), &coils[k], r_mm);
        }
    }
    println!(
        "  every coil watertight and manifold: {all_ok}; worst volume error {worst_v:.2e}, worst surface error {worst_s:.3} mm"
    );

    // Coil-to-coil clearance between centrelines, and the gap left between conductors.
    let samples: Vec<Vec<[f64; 3]>> = coils.iter().map(|cc| cc.sample(512)).collect();
    let mut gap = f64::INFINITY;
    for k in 0..coils.len() {
        let j = (k + 1) % coils.len();
        gap = gap.min(min_distance(&samples[k], &samples[j]));
    }
    println!(
        "closest neighbouring centrelines: {:.1} mm, so {:.1} mm between {:.0} mm conductors",
        gap * MM,
        gap * MM - 2.0 * r_mm,
        2.0 * r_mm
    );

    let all = merge(&parts.iter().map(|p| p.mesh.clone()).collect::<Vec<_>>());
    let path = out.join("coils.stl");
    save_mesh(&all, &path).expect("save coils.stl");
    println!("wrote {} ({} faces)", path.display(), all.faces.len());

    // A look at the set: one field period in copper, the rest grey.
    let per = c.coils_per_period;
    let scene = render::Scene {
        items: parts
            .iter()
            .enumerate()
            .map(|(k, p)| (p.mesh.clone(), if k < per { [205, 120, 60] } else { [170, 170, 175] }))
            .collect(),
    };
    let path = out.join("coils.bmp");
    render::render_bmp(&scene, [-0.55, -0.45, -0.7], 1400, 1000, &path);
    println!("wrote {}", path.display());
}

/// `mesh-probe`: what makes the octree meshes non-manifold, and where their surface error sits.
pub fn probe() {
    let c = design::baseline();
    let coil = c.coil(1);
    let r_mm = design::WIRE_RADIUS * MM;
    let s = Solid::pipe_spline(control_points(&coil, 64, true), r_mm);
    let show = |name: &str, m: &IndexedMesh, secs: f64| {
        let r = validate_mesh(m);
        let line: Vec<[f64; 3]> = coil.sample(8192).into_iter().map(|q| q.map(|x| x * MM)).collect();
        let mut d: Vec<f64> = m
            .vertices
            .par_iter()
            .map(|v| {
                let mut d2 = f64::INFINITY;
                for q in &line {
                    d2 = d2.min((v.x - q[0]).powi(2) + (v.y - q[1]).powi(2) + (v.z - q[2]).powi(2));
                }
                (d2.sqrt() - r_mm).abs()
            })
            .collect();
        d.sort_by(f64::total_cmp);
        let pct = |p: f64| d[((d.len() - 1) as f64 * p) as usize];
        let over = d.iter().filter(|x| **x > 0.1).count();
        println!(
            "  {name:<34} {:>7} faces  boundary {} non-manifold {} degenerate {} duplicate {}  |dev| median {:.4} p99 {:.4} max {:.3} mm, {over} vertices > 0.1 mm  {secs:.1} s",
            m.faces.len(),
            r.boundary_edge_count,
            r.non_manifold_edge_count,
            r.degenerate_face_count,
            r.duplicate_face_count,
            pct(0.5),
            pct(0.99),
            d[d.len() - 1]
        );
    };
    for cell in [1.0, 0.5] {
        let t0 = Instant::now();
        let raw = s.mesh_adaptive_par(cell).geometry;
        let secs = t0.elapsed().as_secs_f64();
        show(&format!("adaptive_par {cell} mm, raw"), &raw, secs);
        let mut w = raw.clone();
        weld_vertices(&mut w, 1e-3);
        w.faces.retain(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
        remove_unreferenced_vertices(&mut w);
        show(&format!("adaptive_par {cell} mm, welded"), &w, secs);
        let t0 = Instant::now();
        let seq = s.mesh_adaptive(cell).geometry;
        show(&format!("adaptive (sequential) {cell} mm"), &seq, t0.elapsed().as_secs_f64());
    }
}

/// `sdf-probe`: is pipe_spline's distance right where the octree mesh is off?
pub fn sdf_probe() {
    let c = design::baseline();
    let coil = c.coil(1);
    let r_mm = design::WIRE_RADIUS * MM;
    let line: Vec<[f64; 3]> = coil.sample(8192).into_iter().map(|q| q.map(|x| x * MM)).collect();
    let true_sdf = |v: &Point3<f64>| {
        let mut d2 = f64::INFINITY;
        for q in &line {
            d2 = d2.min((v.x - q[0]).powi(2) + (v.y - q[1]).powi(2) + (v.z - q[2]).powi(2));
        }
        d2.sqrt() - r_mm
    };
    for n in [64, 256] {
        let s = Solid::pipe_spline(control_points(&coil, n, true), r_mm);
        let m = s.mesh_adaptive_par(0.5).geometry;
        // Vertices more than 0.1 mm off the true tube: what does the solid's own SDF say there?
        let rows: Vec<(f64, f64)> = m
            .vertices
            .par_iter()
            .map(|v| (true_sdf(v), s.evaluate(v)))
            .collect();
        let bad: Vec<&(f64, f64)> = rows.iter().filter(|r| r.0.abs() > 0.1).collect();
        let sdf_err = rows.iter().map(|r| r.1.abs()).fold(0.0, f64::max);
        let bad_sdf = bad.iter().map(|r| r.1.abs()).fold(0.0, f64::max);
        let agree = bad.iter().map(|r| (r.0 - r.1).abs()).fold(0.0, f64::max);
        println!(
            "  {n} control points: {} of {} vertices > 0.1 mm off; |own SDF| at all vertices <= {sdf_err:.4} mm, at the bad ones <= {bad_sdf:.4}; own SDF vs true there differ by <= {agree:.4} mm",
            bad.len(),
            m.vertices.len()
        );
        // Random points near the tube: own SDF against the true distance.
        let mut worst = (0.0_f64, [0.0; 3]);
        let mut seed = 12345_u64;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..20000 {
            let t = rnd() * 2.0 * PI;
            let p = coil.point(t).map(|x| x * MM);
            let q = Point3::new(p[0] + (rnd() - 0.5) * 12.0, p[1] + (rnd() - 0.5) * 12.0, p[2] + (rnd() - 0.5) * 12.0);
            let e = (s.evaluate(&q) - true_sdf(&q)).abs();
            if e > worst.0 {
                worst = (e, [q.x, q.y, q.z]);
            }
        }
        println!("    20,000 random points within 6 mm: own SDF vs true worst {:.4} mm", worst.0);
    }
}
