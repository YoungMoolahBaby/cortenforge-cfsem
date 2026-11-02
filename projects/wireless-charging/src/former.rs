//! `former`: a printable former for the three-coil pad from `layout` (cf-design), in two parts.
//!
//! - `pad_base`: a plate with the two lower coils' channels (open on top), a 1.0 mm floor (a
//!   ferrite sheet goes under it, 1.0 mm from the lower coils, as Qi A6 asks), lead holes.
//! - `pad_lid`: the top coil's channel on its underside, its inner lead's slot running out over
//!   the turns to a hole, and the interface surface on top, 2.2 mm over the top coil.
//!
//! Dowel pins (2 mm) align the two. Units are mm; the assembled frame has the base's top face at
//! z = 0. Each coil's leads go straight down through holes under its inner and outer turns, so
//! only the top coil's inner lead has to cross its own turns (in the slot).

use std::time::Instant;

use cortenforge::cf_design::Solid;
use cortenforge::mesh::io::save_mesh;
use cortenforge::mesh::printability::{PrinterConfig, validate_for_printing};
use cortenforge::mesh::repair::{remove_unreferenced_vertices, validate_mesh, weld_vertices};
use nalgebra::{Point3, UnitQuaternion, Vector3};

use crate::layout::direct_area;

/// Litz diameter (mm) and the channel's clearance on each side.
const WIRE: f64 = 1.15;
const CLEAR: f64 = 0.05;
const FLOOR: f64 = 1.0;
const CHANNEL: f64 = 1.2;
const LID: f64 = 3.4;
const SLOT: (f64, f64) = (1.4, 1.2); // width, height over the channel
const HOLE: f64 = 0.8;
const MARGIN: f64 = 5.0;

#[derive(Clone, Copy)]
struct Pad {
    r_in: f64,
    r_out: f64,
    n: usize,
    s: f64,
}

impl Pad {
    fn pitch(&self) -> f64 {
        (self.r_out - self.r_in) / (self.n - 1) as f64
    }
    /// Channel radii: the inner and outer turns' centres less and plus half a pitch and clearance.
    fn channel(&self) -> (f64, f64) {
        let h = 0.5 * self.pitch().max(WIRE) + CLEAR;
        (self.r_in - h, self.r_out + h)
    }
    fn half(&self) -> (f64, f64) {
        let c = self.channel().1;
        (self.s + c + MARGIN, c + MARGIN)
    }
    fn inner_lead(&self) -> f64 {
        self.r_out + 2.5
    }
    fn dowels(&self) -> [(f64, f64); 2] {
        let (x, y) = self.half();
        [(x - 3.5, y - 3.5), (-(x - 3.5), -(y - 3.5))]
    }
}

fn bx(c: (f64, f64, f64), size: (f64, f64, f64)) -> Solid {
    Solid::cuboid(Vector3::new(0.5 * size.0, 0.5 * size.1, 0.5 * size.2))
        .translate(Vector3::new(c.0, c.1, c.2))
}
fn cyl(x: f64, y: f64, z0: f64, z1: f64, r: f64) -> Solid {
    Solid::cylinder(r, 0.5 * (z1 - z0)).translate(Vector3::new(x, y, 0.5 * (z0 + z1)))
}
fn annulus(x: f64, y: f64, z0: f64, z1: f64, r0: f64, r1: f64) -> Solid {
    cyl(x, y, z0, z1, r1).subtract(cyl(x, y, z0 - 1.0, z1 + 1.0, r0))
}

fn base(p: &Pad) -> Solid {
    let (hx, hy) = p.half();
    let (c0, c1) = p.channel();
    let mut s = bx(
        (0.0, 0.0, -0.5 * (FLOOR + CHANNEL)),
        (2.0 * hx, 2.0 * hy, FLOOR + CHANNEL),
    );
    for sx in [-p.s, p.s] {
        s = s.subtract(annulus(sx, 0.0, -CHANNEL, 1.0, c0, c1));
        // Leads under the outer and inner turns, on the outward side.
        for r in [p.r_out, p.r_in] {
            s = s.subtract(cyl(
                sx + sx.signum() * r,
                0.0,
                -FLOOR - CHANNEL - 1.0,
                1.0,
                HOLE,
            ));
        }
    }
    // Where the two channels nearly meet (|y| < 5.5 mm) the wall between them is a cusp under
    // 1 mm (zero at y = 0); join them there instead (the outer turns pass 0.05 mm apart anyway).
    s = s.subtract(bx(
        (0.0, 0.0, 0.5 * (1.0 - CHANNEL)),
        (1.2, 11.0, CHANNEL + 1.0),
    ));
    // The top coil's leads: its outer turn at -y, its inner lead's slot end at +y.
    for y in [-p.r_out, p.inner_lead()] {
        s = s.subtract(cyl(0.0, y, -FLOOR - CHANNEL - 1.0, 1.0, HOLE));
    }
    for (x, y) in p.dowels() {
        s = s.subtract(cyl(x, y, -FLOOR - CHANNEL - 1.0, 1.0, 1.05));
    }
    s
}

/// The lid in the assembled frame (channel down, interface up).
fn lid(p: &Pad) -> Solid {
    let (hx, hy) = p.half();
    let (c0, c1) = p.channel();
    let mut s = bx((0.0, 0.0, 0.5 * LID), (2.0 * hx, 2.0 * hy, LID));
    s = s.subtract(annulus(0.0, 0.0, -1.0, CHANNEL, c0, c1));
    // The inner lead's slot, from the inner turn out past the coil, open downward at its end.
    let (y0, y1) = (c0, p.inner_lead() + HOLE);
    s = s.subtract(bx(
        (0.0, 0.5 * (y0 + y1), 0.5 * (CHANNEL + SLOT.1 - 1.0)),
        (SLOT.0, y1 - y0, CHANNEL + SLOT.1 + 1.0),
    ));
    for (x, y) in p.dowels() {
        s = s.subtract(cyl(x, y, -1.0, LID + 1.0, 1.05));
    }
    s
}

/// Mesh, weld, check printability as modelled (z up), save; (faces, watertight, cm^3, verdict).
fn part(name: &str, solid: &Solid, cell: f64, out: &std::path::Path) {
    let t0 = Instant::now();
    let mut m = solid.mesh(cell).geometry;
    weld_vertices(&mut m, 1e-3);
    m.faces
        .retain(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
    remove_unreferenced_vertices(&mut m);
    let r = validate_mesh(&m);
    let verdict = match validate_for_printing(&m, &PrinterConfig::fdm_default()) {
        Ok(v) => {
            let mut kinds: Vec<String> = v
                .issues
                .iter()
                .filter(|i| i.is_critical())
                .map(|i| format!("{:?}", i.issue_type))
                .collect();
            for i in v.issues.iter().filter(|i| i.is_critical()) {
                println!(
                    "      critical: {} at {:?}",
                    i.description,
                    i.location.map(|l| (l.x, l.y, l.z))
                );
            }
            kinds.sort();
            kinds.dedup();
            format!(
                "{} (critical {kinds:?})",
                if v.is_printable() {
                    "printable"
                } else {
                    "not printable"
                }
            )
        }
        Err(e) => format!("check failed: {e}"),
    };
    let vol = m.signed_volume();
    println!(
        "   {name:<9} {:>8} faces, watertight {}, {:.1} cm^3 ({:.0} g solid PLA), {verdict}, {:.1} s",
        m.faces.len(),
        r.is_watertight && r.is_manifold,
        vol / 1e3,
        vol / 1e3 * 1.24,
        t0.elapsed().as_secs_f64()
    );
    save_mesh(&m, out.join(format!("{name}.stl"))).expect("write STL");
}

pub fn run() {
    // 1. The layout's optimum, made exact, and a variant with slack between turns.
    let opt = Pad {
        r_in: 6.15,
        r_out: 25.70,
        n: 18,
        s: 26.30,
    };
    // The Qi A11 winds the same Litz at a 1.178 mm pitch.
    let slack = Pad {
        r_in: 25.70 - 17.0 * 1.178,
        r_out: 25.70,
        n: 18,
        s: 26.30,
    };
    println!("1. The three-coil pad (layout's score, exact rings, no table):");
    let mut areas = Vec::new();
    for (name, p) in [("as optimized", opt), ("A11 pitch", slack)] {
        let (hard, soft, l, k0, q1) = direct_area(p.r_in * 1e-3, p.r_out * 1e-3, p.n, p.s * 1e-3);
        println!(
            "   {name:<13} r {:.2}-{:.2} mm, {} turns at {:.2} mm, {:.2} mm apart: {hard:.0} mm^2 (soft {soft:.1}); L {:.2} uH, centred k {k0:.3}, Q1 {q1:.0}",
            p.r_in,
            p.r_out,
            p.n,
            p.pitch(),
            p.s,
            l * 1e6
        );
        areas.push(soft);
    }
    // Chosen: the A11 pitch. At 1.15 mm the 18 turns of 1.15 mm Litz would have 0.1 mm of play
    // between them in all, which can't be wound reliably; Qi's own coil uses 1.178 mm.
    let p = slack;
    let (c0, c1) = p.channel();
    println!(
        "   chosen: the A11 pitch ({:+.1}% of the soft score, for {:.3} mm of play between turns); channels r {c0:.2}-{c1:.2} mm, \
         {:.2} mm wide for {} turns of {WIRE} mm",
        (areas[1] / areas[0] - 1.0) * 100.0,
        p.pitch() - WIRE,
        c1 - c0,
        p.n
    );

    // 2. Parts.
    let (b, l) = (base(&p), lid(&p));
    let probe = |s: &Solid, x: f64, y: f64, z: f64| s.evaluate(&Point3::new(x, y, z));
    let (hx, hy) = p.half();
    println!(
        "2. Former {:.1} x {:.1} mm; base {:.1} mm (floor {FLOOR}), lid {LID} mm (interface {:.2} mm over the top coil)",
        2.0 * hx,
        2.0 * hy,
        FLOOR + CHANNEL,
        LID - WIRE
    );
    println!(
        "   probes (~0 on a surface): base channel floor {:+.3}, inner wall {:+.3}, outer wall {:+.3}; lid channel ceiling {:+.3}, \
         slot ceiling {:+.3}, interface {:+.3}",
        probe(&b, p.s, 15.0, -CHANNEL),
        probe(&b, p.s, c0, -0.5 * CHANNEL),
        probe(&b, p.s, c1, -0.5 * CHANNEL),
        probe(&l, 0.0, -15.0, CHANNEL),
        probe(&l, 0.0, 15.0, CHANNEL + SLOT.1),
        probe(&l, -20.0, -20.0, LID)
    );
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("out")
        .join("parts");
    std::fs::create_dir_all(&out).expect("create out/parts");
    part("pad_base", &b, 0.2, &out);
    // The lid prints upside down: interface on the bed, channel and slot open upward.
    let flip = UnitQuaternion::from_axis_angle(&Vector3::x_axis(), std::f64::consts::PI);
    part("pad_lid", &l.rotate(flip), 0.2, &out);
    println!("   STL in {}", out.display());
}
