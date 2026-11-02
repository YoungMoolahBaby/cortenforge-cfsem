//! `parts`: the demonstrator's printable parts (cf-design), checked for printability and exported,
//! with the landing-skid window from the force table.
//!
//! Units are mm in this module. Track frame: x along the track, z up, the steel rail's underside
//! at z = 0; the stator board (1.5 mm) below it; the cart top `gap` below the steel.

use std::time::Instant;

use cortenforge::cf_design::Solid;
use cortenforge::mesh::io::save_mesh;
use cortenforge::mesh::printability::{PrinterConfig, validate_for_printing};
use cortenforge::mesh::repair::{remove_unreferenced_vertices, validate_mesh, weld_vertices};
use cortenforge::mesh::types::IndexedMesh;
use nalgebra::Vector3;

use crate::corner::{Table, corner_coil};

const G: f64 = 9.81;
const CORNER_KG: f64 = 0.125;
/// Cart body (mm).
const BODY: (f64, f64, f64) = (140.0, 80.0, 30.0);
const CORNERS: [(f64, f64); 4] = [(-50.0, -25.0), (50.0, -25.0), (50.0, 25.0), (-50.0, 25.0)];
/// Magnet + coil assembly pocket: coil outer radius 9.54 mm plus 0.2 mm clearance; 4.0 mm deep.
const POCKET: (f64, f64) = (9.75, 4.0);
const DECK_T: f64 = 6.0;
const WALL: f64 = 2.5;
/// Gap sensor holes: 8 mm, beside each magnet towards the cart's middle (no room outboard: the
/// pocket reaches 59.75 mm and the body ends at 70 mm).
const SENSOR: (f64, f64) = (4.0, 16.0); // radius, x offset inward from the magnet
const BOARD: f64 = 1.5; // stator board thickness under the steel (mm)
const STEEL: (f64, f64) = (60.0, 3.0); // flat bar width, thickness
const CLEAR: f64 = 2.0; // side clearance to the guide lips
const REST_GAP: f64 = 6.5; // mm, where the station shelf holds the cart (lift-off start)

fn bx(cx: f64, cy: f64, cz: f64, sx: f64, sy: f64, sz: f64) -> Solid {
    Solid::cuboid(Vector3::new(0.5 * sx, 0.5 * sy, 0.5 * sz)).translate(Vector3::new(cx, cy, cz))
}

fn cyl(cx: f64, cy: f64, z0: f64, z1: f64, r: f64) -> Solid {
    Solid::cylinder(r, 0.5 * (z1 - z0)).translate(Vector3::new(cx, cy, 0.5 * (z0 + z1)))
}

/// The deck, printed flat (bottom on the bed, pockets and skids up): a plate with the four magnet
/// and coil pockets, the gap-sensor holes, and the landing skids.
fn deck(skid: f64) -> Solid {
    let mut d = bx(0.0, 0.0, 0.5 * DECK_T, BODY.0, BODY.1, DECK_T);
    for &(x, y) in &CORNERS {
        d = d.subtract(cyl(x, y, DECK_T - POCKET.1, DECK_T + 1.0, POCKET.0));
        let sx = x - x.signum() * SENSOR.1;
        d = d.subtract(cyl(sx, y, -1.0, DECK_T + 1.0, SENSOR.0));
        // Wire passage from the pocket floor down through the deck.
        d = d.subtract(cyl(x, y, -1.0, DECK_T, 1.5));
    }
    // Four skids near the ends, inside the stator board's width.
    for (x, y) in [(-62.0, -12.0), (-62.0, 12.0), (62.0, -12.0), (62.0, 12.0)] {
        d = d.union(cyl(x, y, DECK_T - 0.5, DECK_T + skid, 3.0));
    }
    d
}

/// The hull: an open tray for battery and electronics, printed open side up.
fn hull() -> Solid {
    let h = BODY.2 - DECK_T;
    bx(0.0, 0.0, 0.5 * h, BODY.0, BODY.1, h).subtract(bx(
        0.0,
        0.0,
        0.5 * h + WALL,
        BODY.0 - 2.0 * WALL,
        BODY.1 - 2.0 * WALL,
        h,
    ))
}

/// The rail hanger's cross-section (y, z), extruded `len` along x and printed standing on end:
/// a top plate over the steel bar, tabs under the stator board's edges, side lips guiding the
/// cart, and at a station, ledges where the cart rests at REST_GAP.
fn hanger(len: f64, station: bool) -> Solid {
    let lip_in = 0.5 * BODY.1 + CLEAR; // inner face of the lips
    let lip_t = 3.0;
    let top = (STEEL.1, STEEL.1 + 5.0);
    let bottom = -(REST_GAP + BODY.2) - if station { 3.0 } else { -6.0 };
    let slab = |y0: f64, y1: f64, z0: f64, z1: f64| {
        bx(0.0, 0.5 * (y0 + y1), 0.5 * (z0 + z1), len, y1 - y0, z1 - z0)
    };
    let mut s = slab(-lip_in - lip_t, lip_in + lip_t, top.0, top.1);
    for sign in [-1.0, 1.0] {
        let (a, b) = (sign * lip_in, sign * (lip_in + lip_t));
        s = s.union(slab(a.min(b), a.max(b), bottom, top.1));
        // Tabs holding the steel bar and board: from the board edge out to the lip.
        let (c, d) = (sign * (0.5 * STEEL.0 - 1.0), sign * lip_in);
        s = s.union(slab(c.min(d), c.max(d), -BOARD - 1.5, -BOARD));
        // Spacer beside the steel bar.
        let (e, f) = (sign * 0.5 * STEEL.0, sign * lip_in);
        s = s.union(slab(e.min(f), e.max(f), -BOARD, top.0));
        if station {
            // Ledge under the cart's edge: its top is where the cart's bottom rests.
            let (g, h) = (sign * (lip_in - 6.0), sign * lip_in);
            let z_rest = -(REST_GAP + BODY.2);
            s = s.union(slab(g.min(h), g.max(h), bottom, z_rest));
        }
    }
    // Two mounting holes through the top plate, diamond-shaped (a square turned 45 deg about the
    // hole axis): printed on end the hole runs sideways, and a round hole's top would be a flat
    // overhang; a diamond's faces are all at 45 deg.
    let turn =
        nalgebra::UnitQuaternion::from_axis_angle(&Vector3::z_axis(), 0.25 * std::f64::consts::PI);
    for x in [-0.25 * len, 0.25 * len] {
        let hole = Solid::cuboid(Vector3::new(2.2, 2.2, 0.5 * (top.1 - top.0) + 1.0))
            .rotate(turn)
            .translate(Vector3::new(x, 0.0, 0.5 * (top.0 + top.1)));
        s = s.subtract(hole);
    }
    // Printed standing on end: turn the extrusion axis (x) to z.
    s.rotate(nalgebra::UnitQuaternion::from_axis_angle(
        &Vector3::y_axis(),
        -0.5 * std::f64::consts::PI,
    ))
}

/// Mesh, weld, check printability as modelled (z up), save; returns the welded mesh.
fn part(name: &str, solid: &Solid, cell: f64, out: &std::path::Path) -> IndexedMesh {
    let t0 = Instant::now();
    let mut m = solid.mesh(cell).geometry;
    weld_vertices(&mut m, 1e-3);
    m.faces
        .retain(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
    remove_unreferenced_vertices(&mut m);
    let r = validate_mesh(&m);
    let vol = m.signed_volume();
    let verdict = match validate_for_printing(&m, &PrinterConfig::fdm_default()) {
        Ok(v) => {
            let crit: Vec<String> = v
                .issues
                .iter()
                .filter(|i| i.is_critical())
                .map(|i| format!("{:?}", i.issue_type))
                .collect();
            let mut kinds = crit.clone();
            kinds.sort();
            kinds.dedup();
            let area: f64 = v.overhangs.iter().map(|o| o.area).sum();
            let (zmin, _) = m
                .vertices
                .iter()
                .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.z), a.1.max(p.z)));
            for i in v.issues.iter().filter(|i| i.is_critical()) {
                if let Some(l) = i.location {
                    println!(
                        "      critical at ({:.1}, {:.1}, {:.1}) mm, {:.1} above the bed: {}",
                        l.x,
                        l.y,
                        l.z,
                        l.z - zmin,
                        i.description
                    );
                }
            }
            format!(
                "{} (critical: {} {:?}; overhang area {:.1} mm^2)",
                if v.is_printable() {
                    "printable"
                } else {
                    "not printable"
                },
                crit.len(),
                kinds,
                area
            )
        }
        Err(e) => format!("check failed: {e}"),
    };
    println!(
        "   {name:<16} {:>8} faces, watertight {}, {:>7.1} cm^3 ({:>5.0} g solid PLA), {verdict}, {:.1} s",
        m.faces.len(),
        r.is_watertight && r.is_manifold,
        vol / 1e3,
        vol / 1e3 * 1.24,
        t0.elapsed().as_secs_f64()
    );
    save_mesh(&m, out.join(format!("{name}.stl"))).expect("write STL");
    m
}

pub fn run() {
    // 1. The landing-skid window from the force table.
    let tab = Table::build(&corner_coil());
    let w = CORNER_KG * G;
    let push_off = |g: f64| (tab.eval(&tab.a, g).0 - w) / tab.eval(&tab.b, g).0;
    println!(
        "1. Landing skids: on power loss the magnets pull the cart up until the skids touch the stator board."
    );
    println!(
        "   Current to push a corner back off (unloaded, the hardest case) at the landed gap:"
    );
    for g in [3.5, 4.0, 4.5, 5.0, 5.5] {
        println!(
            "     landed gap {g:.1} mm: {:.2} A (limit 10 A)",
            push_off(g * 1e-3).abs()
        );
    }
    // Landed gap where push-off needs 5 A (half the limit), for the record.
    let (mut lo, mut hi) = (0.003, 0.0059);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if push_off(mid).abs() > 5.0 {
            lo = mid
        } else {
            hi = mid
        }
    }
    let g_half = hi * 1e3;
    let min_run = 5.034; // mm, the lowest corner gap in `cart` (front-left, 100 g payload)
    println!(
        "   push-off within 5 A needs a landed gap >= {g_half:.2} mm, but running gaps go down to {min_run:.3} mm \
         (`cart`, 100 g off-centre): no room at that margin"
    );
    // Chosen: land 0.5 mm below the lowest running gap, and accept the push-off current there.
    let g_land = min_run - 0.5;
    let skid = g_land - BOARD;
    println!(
        "   chosen: landed gap {g_land:.2} mm (0.5 mm under the lowest running gap): push-off {:.2} A of the 10 A \
         limit; skids {skid:.2} mm tall (the board is {BOARD} mm)",
        push_off(g_land * 1e-3).abs()
    );

    // 2. Parts.
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("out")
        .join("parts");
    std::fs::create_dir_all(&out).expect("create out/parts");
    println!(
        "2. Printable parts (cf-design, marching cubes, welded; FDM default check as modelled, z up):"
    );
    let d = deck(skid);
    let probe = |s: &Solid, x: f64, y: f64, z: f64| s.evaluate(&nalgebra::Point3::new(x, y, z));
    println!(
        "   deck probes (mm, ~0 on a surface): pocket floor {:+.3}, pocket wall {:+.3}, skid top {:+.3}",
        probe(&d, 55.0, 25.0, DECK_T - POCKET.1),
        probe(&d, 50.0 + POCKET.0, 25.0, DECK_T - 2.0),
        probe(&d, 62.0, 12.0, DECK_T + skid)
    );
    let deck_m = part("cart_deck", &d, 0.25, &out);
    let hull_m = part("cart_hull", &hull(), 0.4, &out);
    part("rail_hanger", &hanger(100.0, false), 0.4, &out);
    part("station_hanger", &hanger(100.0, true), 0.4, &out);
    println!("   STL in {}", out.display());
    let printed_g = (deck_m.signed_volume() + hull_m.signed_volume()) / 1e3 * 1.24;
    println!(
        "3. Mass budget (cart, 500 g in the model with its 200 g nominal payload): printed deck + hull {printed_g:.0} g \
         solid PLA (about {:.0} g at 40% effective density); magnets 4 x {:.1} g; coils 4 x {:.1} g copper",
        printed_g * 0.4,
        7.5e-3 * std::f64::consts::PI * 0.6 * 0.6 * 0.4 * 1e3,
        corner_coil().turns() as f64
            * 2.0
            * std::f64::consts::PI
            * 0.8
            * 0.0321
            * 0.0321
            * std::f64::consts::PI
            / 4.0
            * 8.96
    );
}
