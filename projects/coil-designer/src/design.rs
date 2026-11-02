//! `design`: the least-power Helmholtz-type pair that meets a robust uniformity spec, and its
//! printed coil former.
//!
//! Spec: a 10 mm diameter sphere uniform to 1e-4 at 1 mT, holding for any spacing error within
//! +-0.1 mm (FDM print tolerance), at least 5 mm between the packs. Each coil is exactly 100 turns
//! of AWG 22 at 50 mm mean radius, layer-wound with a partial last layer, so power compares
//! shapes at equal copper. For each axial turn count the search finds the spacing window that
//! meets the spec, centres the design in it, and keeps the window only if it is at least
//! +-0.1 mm wide; the least power at 1 mT wins.
//!
//! Step 2 (`optimize`) found a grid over the integer with a 1-D search inside to be the right
//! tool for this problem shape, so this uses that, not sim-opt.

use std::f64::consts::PI;
use std::time::Instant;

use cortenforge::cf_design::Solid;
use cortenforge::mesh::io::save_stl;
use cortenforge::mesh::printability::{
    PrinterConfig, find_optimal_orientation, place_on_build_plate, validate_for_printing,
};
use cortenforge::mesh::repair::{remove_unreferenced_vertices, validate_mesh, weld_vertices};
use nalgebra::Vector3;
use rayon::prelude::*;

use crate::coil::{Coil, Pack, Wire, uniform_radius};

const R_MEAN: f64 = 0.05; // m
const TURNS: usize = 100; // per coil, exactly
const EPS: f64 = 1e-4;
const RHO_REQ: f64 = 5e-3; // m, radius of the sphere that must stay uniform
const TOL: f64 = 0.1e-3; // m, spacing error the design must tolerate
const B_TARGET: f64 = 1e-3; // T
const GAP_MIN: f64 = 5e-3; // m, two 2.5 mm former walls
const PER_LAYER: (usize, usize) = (4, 100); // axial turns per layer

// Former, in mm.
const WALL: f64 = 2.5; // tube wall and flange thickness
const CLEARANCE: f64 = 0.2; // added to each channel's width for winding slack
const RIM: f64 = 1.5; // flange height above the full winding

/// One coil's pack: exactly TURNS turns, `per_layer` axially, last layer partial.
fn pack(per_layer: usize, z_mid: f64) -> Pack {
    let wire = Wire::awg(22);
    let layers = TURNS.div_ceil(per_layer);
    let last = TURNS - (layers - 1) * per_layer;
    let depth = layers as f64 * wire.outer;
    Pack {
        r_in: R_MEAN - 0.5 * depth,
        z_mid,
        layers,
        per_layer,
        wire,
        sign: 1.0,
        partial: if last == per_layer { 0 } else { last },
    }
}

fn pair(per_layer: usize, spacing: f64) -> Coil {
    Coil {
        packs: vec![
            pack(per_layer, -0.5 * spacing),
            pack(per_layer, 0.5 * spacing),
        ],
    }
}

/// (m) EPS-uniform sphere radius.
fn rho(per_layer: usize, spacing: f64) -> f64 {
    let c = pair(per_layer, spacing);
    let f = |r: &[f64], z: &[f64]| c.field(1.0, r, z);
    uniform_radius(&f, EPS, 0.3 * R_MEAN)
}

/// A shape that meets the spec somewhere: its spacing window and centred design.
struct Candidate {
    per_layer: usize,
    lo: f64,
    hi: f64,
    spacing: f64,
    rho_mid: f64,
    rho_edges: f64,
    power: f64,
    evals: usize,
}

/// Spacing window where rho >= RHO_REQ, around the best point of a coarse scan.
fn window(per_layer: usize) -> Option<Candidate> {
    let width = pack(per_layer, 0.0).width();
    let s_min = (width + GAP_MIN + TOL).max(0.9 * R_MEAN);
    let s_max = 1.6 * R_MEAN;
    let n = ((s_max - s_min) / (0.005 * R_MEAN)).ceil() as usize;
    let mut evals = 0;
    let grid: Vec<(f64, f64)> = (0..=n)
        .map(|k| {
            let s = s_min + (s_max - s_min) * k as f64 / n as f64;
            evals += 1;
            (s, rho(per_layer, s))
        })
        .collect();
    let (k_best, &(s_best, r_best)) = grid
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.1.total_cmp(&b.1.1))?;
    // Refine the peak: the ridge is often narrower than the coarse step.
    let mut peak = (s_best, r_best);
    let step = (s_max - s_min) / n as f64;
    for k in -50..=50 {
        let s = (s_best + step * k as f64 / 50.0).clamp(s_min, s_max);
        let r = rho(per_layer, s);
        evals += 1;
        if r > peak.1 {
            peak = (s, r);
        }
    }
    if peak.1 < RHO_REQ {
        return None;
    }
    // Bisect each side for the crossing, from the peak out to the first grid point below.
    let mut edge = |dir: f64| -> f64 {
        let mut inside = peak.0;
        let mut k = k_best as isize;
        let outside = loop {
            k += dir as isize;
            if k < 0 || k as usize > n {
                return if dir < 0.0 { s_min } else { s_max };
            }
            let (s, r) = grid[k as usize];
            if r < RHO_REQ && (s - peak.0) * dir > 0.0 {
                break s;
            }
        };
        let mut out = outside;
        for _ in 0..30 {
            let mid = 0.5 * (inside + out);
            evals += 1;
            if rho(per_layer, mid) >= RHO_REQ {
                inside = mid
            } else {
                out = mid
            }
        }
        inside
    };
    let lo = edge(-1.0);
    let hi = edge(1.0);
    let spacing = 0.5 * (lo + hi);
    let c = pair(per_layer, spacing);
    let b0 = c.field(1.0, &[0.0], &[0.0]).1[0];
    let power = (B_TARGET / b0).powi(2) * c.resistance();
    let rho_edges = rho(per_layer, spacing - TOL).min(rho(per_layer, spacing + TOL));
    evals += 3;
    Some(Candidate {
        per_layer,
        lo,
        hi,
        spacing,
        rho_mid: rho(per_layer, spacing),
        rho_edges,
        power,
        evals,
    })
}

/// Former dimensions in mm, for one half printed standing on its mid-plane face (z = 0, the
/// bed). The two halves are identical; they butt at the mid-plane, aligned by dowels.
struct Dims {
    r_bore: f64,
    r_in: f64,
    r_flange: f64,
    /// Top of the winding floor: the wound pack rests here, so it sets the spacing.
    floor: f64,
    /// Top of the channel (pack width + clearance above the floor).
    channel_top: f64,
    /// Top of the tube: one wall above the channel, for the washer to slide on.
    top: f64,
}

/// Collar at the bed, radial reach past the winding surface (mm). Its height puts its outer top
/// edge on the support cone, so the cone runs straight into the collar wall with no crevice.
const COLLAR: f64 = 4.0;

/// Dowel holes in the collar: count, radius, depth before the pointed tip (mm). Sized for
/// 1.75 mm filament as dowels.
const DOWEL: (usize, f64, f64) = (3, 0.95, 2.5);
/// Rise per unit radius of the support cone under the floor flange: 1.2 is 40 deg from vertical.
const CONE_RISE: f64 = 1.2;
/// Washer bore clearance over the tube (mm).
const WASHER_FIT: f64 = 0.15;
/// Print layer height (mm). The winding floor is snapped to a layer boundary: a slicer can only
/// put a horizontal face on one, and the floor sets the spacing.
const LAYER_H: f64 = 0.1;

fn dims(per_layer: usize, spacing: f64) -> Dims {
    let p = pack(per_layer, 0.0);
    let r_in = p.r_in * 1e3;
    let w = p.width() * 1e3;
    let floor = ((0.5 * spacing * 1e3 - 0.5 * w) / LAYER_H).round() * LAYER_H;
    let channel_top = floor + w + CLEARANCE;
    Dims {
        r_bore: r_in - WALL,
        r_in,
        r_flange: r_in + p.depth() * 1e3 + RIM,
        floor,
        channel_top,
        top: channel_top + WALL,
    }
}

/// A z-aligned annulus between radii, spanning z0..z1.
fn ring(r0: f64, r1: f64, z0: f64, z1: f64) -> Solid {
    let h = 0.5 * (z1 - z0);
    Solid::cylinder(r1, h)
        .subtract(Solid::cylinder(r0, h + 1.0))
        .translate(Vector3::new(0.0, 0.0, z0 + h))
}

/// One half: tube, floor flange with a support cone under it, bed collar with dowel holes.
fn former_half(d: &Dims) -> Solid {
    let flange_bottom = d.floor - WALL;
    let cone_bottom = flange_bottom - CONE_RISE * (d.r_flange - d.r_in);
    assert!(cone_bottom > 0.0, "no room for the support cone in the gap");
    // Cone with its apex down and radius r_flange at flange_bottom, clipped to its slab.
    let apex = flange_bottom - CONE_RISE * d.r_flange;
    let flip = nalgebra::UnitQuaternion::from_axis_angle(&Vector3::x_axis(), PI);
    let cone = Solid::cone(d.r_flange, CONE_RISE * d.r_flange)
        .rotate(flip)
        .translate(Vector3::new(0.0, 0.0, apex))
        .intersect(
            Solid::cylinder(d.r_flange + 1.0, 0.5 * (flange_bottom - cone_bottom))
                .translate(Vector3::new(0.0, 0.0, 0.5 * (flange_bottom + cone_bottom))),
        );
    let mut body = ring(d.r_bore, d.r_in, 0.0, d.top)
        .union(ring(d.r_bore, d.r_flange, flange_bottom, d.floor))
        .union(cone)
        .union(ring(
            d.r_bore,
            d.r_in + COLLAR,
            0.0,
            cone_bottom + CONE_RISE * COLLAR,
        ))
        .subtract(Solid::cylinder(d.r_bore, d.top + 2.0));
    // Blind dowel holes from the bed, each with a pointed tip so its end needs no bridge.
    let r_hole = 0.5 * (d.r_bore + d.r_in + COLLAR);
    for k in 0..DOWEL.0 {
        let a = 2.0 * PI * k as f64 / DOWEL.0 as f64;
        let at = Vector3::new(r_hole * a.cos(), r_hole * a.sin(), 0.0);
        let shaft = Solid::cylinder(DOWEL.1, 0.5 * DOWEL.2 + 0.5)
            .translate(at + Vector3::new(0.0, 0.0, 0.5 * DOWEL.2 - 0.5));
        let tip = Solid::cone(DOWEL.1, CONE_RISE * DOWEL.1 * 1.5)
            .translate(at + Vector3::new(0.0, 0.0, DOWEL.2 + CONE_RISE * DOWEL.1 * 1.5));
        body = body.subtract(shaft.union(tip));
    }
    body
}

/// The washer that closes a channel: printed flat, slid over the tube end onto the winding.
fn former_washer(d: &Dims) -> Solid {
    ring(d.r_in + WASHER_FIT, d.r_flange, 0.0, WALL)
}

pub fn run() {
    println!(
        "Spec: {:.0} mm sphere uniform to {EPS:.0e} at {:.0} mT, robust to +-{:.1} mm spacing, \
         gap >= {:.0} mm; {TURNS} turns AWG 22 per coil at {:.0} mm mean radius",
        2.0 * RHO_REQ * 1e3,
        B_TARGET * 1e3,
        TOL * 1e3,
        GAP_MIN * 1e3,
        R_MEAN * 1e3
    );
    let t0 = Instant::now();
    let cands: Vec<Option<Candidate>> = (PER_LAYER.0..=PER_LAYER.1)
        .into_par_iter()
        .map(window)
        .collect();
    let evals: usize = cands.iter().flatten().map(|c| c.evals).sum();
    println!(
        "Searched {} shapes in {:.1} s ({} evaluations for the feasible ones)",
        PER_LAYER.1 - PER_LAYER.0 + 1,
        t0.elapsed().as_secs_f64(),
        evals
    );
    let ok: Vec<&Candidate> = cands
        .iter()
        .flatten()
        .filter(|c| c.hi - c.lo >= 2.0 * TOL && c.rho_edges >= RHO_REQ)
        .collect();
    println!(
        "  meets the uniformity spec somewhere: {}",
        cands.iter().flatten().count()
    );
    println!(
        "  and with a window of at least +-{:.1} mm: {}",
        TOL * 1e3,
        ok.len()
    );
    println!("  turns/layer  layers  window (mm)        width   rho mid  rho +-tol  P(1 mT)");
    for c in cands.iter().flatten() {
        let p = pack(c.per_layer, 0.0);
        let mark = if ok.iter().any(|o| o.per_layer == c.per_layer) {
            ' '
        } else {
            'x'
        };
        println!(
            "{mark} {:>10}  {:>6}  {:6.2} - {:6.2}  {:5.2}   {:5.2}    {:5.2}     {:.3} W",
            c.per_layer,
            p.layers,
            c.lo * 1e3,
            c.hi * 1e3,
            (c.hi - c.lo) * 1e3,
            c.rho_mid * 1e3,
            c.rho_edges * 1e3,
            c.power
        );
    }
    let Some(best) = ok.iter().min_by(|a, b| a.power.total_cmp(&b.power)) else {
        println!("No shape meets the spec.");
        return;
    };

    let p = pack(best.per_layer, 0.0);
    let coil = pair(best.per_layer, best.spacing);
    let b0 = coil.field(1.0, &[0.0], &[0.0]).1[0];
    let i = B_TARGET / b0;
    let r = coil.resistance();
    let l = coil.inductance();
    println!("Chosen: least power among those meeting the spec");
    println!(
        "  each coil: {} turns, {} per layer, {} layers (last {}), {:.2} x {:.2} mm, r_in {:.2} mm",
        p.turns(),
        p.per_layer,
        p.layers,
        if p.partial == 0 {
            p.per_layer
        } else {
            p.partial
        },
        p.depth() * 1e3,
        p.width() * 1e3,
        p.r_in * 1e3
    );
    println!(
        "  spacing {:.3} mm between pack centres ({:.4} R), window {:.2}-{:.2} mm",
        best.spacing * 1e3,
        best.spacing / R_MEAN,
        best.lo * 1e3,
        best.hi * 1e3
    );
    println!(
        "  uniform sphere {:.2} mm radius at the centre spacing, >= {:.2} mm across +-{:.1} mm",
        best.rho_mid * 1e3,
        best.rho_edges * 1e3,
        TOL * 1e3
    );
    println!(
        "  {:.3} mT/A; 1 mT at {:.3} A, {:.2} V, {:.3} W; R {:.3} ohm, L {:.3} mH, {:.1} m of wire",
        b0 * 1e3,
        i,
        i * r,
        best.power,
        r,
        l * 1e3,
        coil.wire_length()
    );

    println!(
        "  former: cargo run --release -p coil-designer -- former {} {:.3}",
        best.per_layer,
        best.spacing * 1e3
    );
    build_former(best.per_layer, best.spacing, 0.25, "mc");
}
/// `former <turns per layer> <spacing mm> [mesh cell mm] [mc|dc|adaptive]`.
pub fn run_former(args: &[String]) {
    let per_layer: usize = args.first().and_then(|a| a.parse().ok()).unwrap_or(7);
    let spacing_mm: f64 = args.get(1).and_then(|a| a.parse().ok()).unwrap_or(49.542);
    let cell: f64 = args.get(2).and_then(|a| a.parse().ok()).unwrap_or(0.25);
    let mesher = args.get(3).map_or("mc", String::as_str);
    build_former(per_layer, spacing_mm * 1e-3, cell, mesher);
}

/// The former: two identical halves and two washers, each cf-design solid -> mesh -> accuracy ->
/// printability -> STL.
fn build_former(per_layer: usize, spacing: f64, cell: f64, mesher: &str) {
    let p = pack(per_layer, 0.0);
    let d = dims(per_layer, spacing);
    let w = p.width() * 1e3;
    println!(
        "Former for {} x {} turns at {:.3} mm spacing (mm): bore {:.2}, winding surface {:.2}, \
         flanges {:.2} in diameter",
        p.per_layer,
        p.layers,
        spacing * 1e3,
        2.0 * d.r_bore,
        2.0 * d.r_in,
        2.0 * d.r_flange
    );
    println!(
        "  each half, from its mid-plane face: floor at {:.3}, channel to {:.3} ({:.2} wide for a \
         {w:.2} pack), tube to {:.3}",
        d.floor,
        d.channel_top,
        d.channel_top - d.floor,
        d.top
    );
    let realized = 2.0 * (d.floor + 0.5 * w);
    println!(
        "  floor snapped to {LAYER_H} mm print layers: spacing from the parts, 2 x (floor + \
         pack/2) = {realized:.3} mm ({:+.3} mm from the design)",
        realized - spacing * 1e3
    );
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out");
    std::fs::create_dir_all(&out).expect("create out/");
    let _ = std::fs::remove_file(out.join("former_oriented.stl"));
    let half = former_half(&d);
    // Probes of the key faces: ~0 on a surface.
    let probe = |x: f64, z: f64| half.evaluate(&nalgebra::Point3::new(x, 0.0, z));
    println!(
        "  SDF probes: floor face {:+.4}, winding surface {:+.4}, bore {:+.4}, mid-plane face {:+.4}",
        probe(d.r_in + 2.0, d.floor),
        probe(d.r_in, 0.5 * (d.floor + d.channel_top)),
        probe(d.r_bore, 0.5 * (d.floor + d.channel_top)),
        probe(d.r_bore + 1.0, 0.0)
    );
    for (name, solid, count) in [("half", half, 2), ("washer", former_washer(&d), 2)] {
        part(&out, name, &solid, count, cell, mesher);
    }
}

/// Mesh one part, measure how far the mesh strays from the solid, check printability as
/// modelled (z up, bed at z = 0), run the orientation search, write the STL.
fn part(out: &std::path::Path, name: &str, solid: &Solid, count: usize, cell: f64, mesher: &str) {
    let t0 = Instant::now();
    let mesh = match mesher {
        "dc" => solid.mesh_dc(cell),
        "adaptive" => solid.mesh_adaptive(cell),
        _ => solid.mesh(cell),
    };
    // The mesher leaves coincident and sub-micron-apart vertices unmerged (slivers where the
    // surface grazes grid corners): weld within 1 um, drop faces that collapse.
    let mut welded = mesh.geometry.clone();
    let before = validate_mesh(&welded);
    let merged = weld_vertices(&mut welded, 1e-3);
    let n_faces = welded.faces.len();
    welded
        .faces
        .retain(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
    let collapsed = n_faces - welded.faces.len();
    remove_unreferenced_vertices(&mut welded);
    let after = validate_mesh(&welded);
    println!(
        "  mesh as emitted: watertight {}, manifold {}, {} boundary edges, {} degenerate faces; \
         welded {merged} vertices, dropped {collapsed} collapsed faces -> watertight {}, \
         manifold {}, {} boundary edges, {} degenerate faces",
        before.is_watertight,
        before.is_manifold,
        before.boundary_edge_count,
        before.degenerate_face_count,
        after.is_watertight,
        after.is_manifold,
        after.boundary_edge_count,
        after.degenerate_face_count
    );
    let geom = &welded;
    let errs: Vec<f64> = geom
        .vertices
        .par_iter()
        .map(|v| solid.evaluate(v).abs())
        .collect();
    let worst = errs.iter().copied().fold(0.0, f64::max);
    let far = errs.iter().filter(|&&e| e > 0.25 * cell).count();
    println!(
        "{name} (print {count}): {mesher} mesh at {cell} mm in {:.1} s, {} faces; worst vertex \
         {worst:.3} mm off the surface, {far} vertices beyond a quarter cell",
        t0.elapsed().as_secs_f64(),
        geom.faces.len()
    );
    let config = PrinterConfig::fdm_default();
    let t0 = Instant::now();
    match validate_for_printing(geom, &config) {
        Ok(v) => {
            println!(
                "  printability (FDM default, as modelled), {:.1} s: {}",
                t0.elapsed().as_secs_f64(),
                v.summary()
            );
            let mut kinds: std::collections::BTreeMap<String, (usize, String)> = Default::default();
            for i in &v.issues {
                let e = kinds
                    .entry(format!("{:?} / {}", i.issue_type, i.severity.as_str()))
                    .or_insert((0, i.description.clone()));
                e.0 += 1;
            }
            for (k, (n, first)) in &kinds {
                println!("    {n:>5} x {k}: e.g. \"{first}\"");
            }
            // Do the "self-intersecting" pairs share a vertex (touching, not crossing)?
            let shared = v
                .self_intersecting
                .iter()
                .filter(|s| {
                    let a = geom.faces[s.face_a as usize];
                    let b = geom.faces[s.face_b as usize];
                    a.iter().any(|i| b.contains(i))
                })
                .count();
            if !v.self_intersecting.is_empty() {
                let s = &v.self_intersecting[0];
                let l = s.approximate_location;
                let ta = geom.faces[s.face_a as usize].map(|i| geom.vertices[i as usize]);
                let tb = geom.faces[s.face_b as usize].map(|i| geom.vertices[i as usize]);
                let within = (0..3)
                    .flat_map(|i| (i + 1..3).map(move |j| (i, j)))
                    .map(|(i, j)| (tb[i] - tb[j]).norm())
                    .fold(f64::INFINITY, f64::min);
                let across = ta
                    .iter()
                    .flat_map(|p| tb.iter().map(move |q| (p - q).norm()))
                    .fold(f64::INFINITY, f64::min);
                println!(
                    "    first pair: shortest edge of face b {within:.2e} mm, closest vertices \
                     across the pair {across:.2e} mm"
                );
                for f in [s.face_a, s.face_b] {
                    let tri = geom.faces[f as usize].map(|i| geom.vertices[i as usize]);
                    let n = (tri[1] - tri[0]).cross(&(tri[2] - tri[0]));
                    println!(
                        "    face {f}: {:?}, normal {:.3?}, SDF at centroid {:+.4}",
                        tri.map(|p| [
                            (p.x * 1e3).round() / 1e3,
                            (p.y * 1e3).round() / 1e3,
                            (p.z * 1e3).round() / 1e3
                        ]),
                        n.normalize().as_slice(),
                        solid.evaluate(&nalgebra::Point3::from(
                            (tri[0].coords + tri[1].coords + tri[2].coords) / 3.0
                        ))
                    );
                }
                // Face areas: are the flagged faces zero-area slivers?
                let area = |f: u32| {
                    let [a, b, c] = geom.faces[f as usize];
                    let (a, b, c) = (
                        geom.vertices[a as usize],
                        geom.vertices[b as usize],
                        geom.vertices[c as usize],
                    );
                    0.5 * (b - a).cross(&(c - a)).norm()
                };
                let tiny = v
                    .self_intersecting
                    .iter()
                    .filter(|s| area(s.face_a).min(area(s.face_b)) < 1e-9)
                    .count();
                let degenerate = (0..geom.faces.len() as u32)
                    .filter(|&f| area(f) < 1e-9)
                    .count();
                // Coincident positions under different indices: touching, not crossing.
                let touching = v
                    .self_intersecting
                    .iter()
                    .filter(|s| {
                        let a = geom.faces[s.face_a as usize];
                        let b = geom.faces[s.face_b as usize];
                        a.iter().any(|&i| {
                            b.iter().any(|&j| {
                                (geom.vertices[i as usize] - geom.vertices[j as usize]).norm()
                                    < 1e-9
                            })
                        })
                    })
                    .count();
                println!(
                    "    self-intersecting pairs reported: {}, {shared} sharing a vertex index, \
                     {touching} with coincident vertex positions, {tiny} involving a zero-area \
                     face; first at r {:.2}, z {:.2}; zero-area faces in the mesh: {degenerate}",
                    v.self_intersecting.len(),
                    l.x.hypot(l.y),
                    l.z
                );
            }
            // Where the critical overhangs are.
            for i in v.issues.iter().filter(|i| i.is_critical()).take(12) {
                if let Some(l) = i.location {
                    println!(
                        "    critical {:?} at r {:.2}, z {:.2}: {}",
                        i.issue_type,
                        l.x.hypot(l.y),
                        l.z,
                        i.description
                    );
                }
            }
        }
        Err(e) => println!("  printability check failed: {e}"),
    }
    let orient = find_optimal_orientation(geom, &config, 24);
    let (axis, angle) = orient
        .rotation
        .axis_angle()
        .map_or((Vector3::z(), 0.0), |(a, t)| (a.into_inner(), t));
    println!(
        "  orientation search: rotate {:.0} deg about ({:.2}, {:.2}, {:.2}), overhang {:.0} mm^2",
        angle * 180.0 / PI,
        axis.x,
        axis.y,
        axis.z,
        orient.overhang_area
    );
    let path = out.join(format!("former_{name}.stl"));
    save_stl(&place_on_build_plate(geom), &path, true).expect("write STL");
    let size = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
    println!("  wrote {} ({:.1} MB)", path.display(), size as f64 / 1e6);
}
