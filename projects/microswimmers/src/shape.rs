//! `shape`: the swimmer's helix (radius, pitch angle, wire, turns) chosen for the fastest climb
//! at 2 mT, printable in resin. Then the winner as a solid (cf-design), meshed, checked
//! (mesh-printability, SLA) and written as STL.
//!
//! The head is set by the magnet and the printer: CortenForge's SLA preset wants walls of at
//! least 0.4 mm, so the 2 x 2 mm magnet sits in a 2.9 mm head, and it goes in through an open,
//! flat nose (sealed after). Climb = top speed at 2 mT (at step-out) less sinking, axis upright.
//!
//! The search is CortenForge's RicherSA (sim-opt) on a dummy one-step model, as in coil-designer:
//! three zero-gear motors carry the design, the reward builds the swimmer and scores it. Every
//! evaluation is logged and the best evaluated design is kept. One run for each whole number of
//! turns. Each evaluation uses the model at 6 and 8 points round the wire, extrapolated; the
//! baseline and the winners are checked at 12 and 16.

use std::f64::consts::PI;
use std::fmt::Write as _;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use cortenforge::cf_design::Solid;
use cortenforge::mesh::io::save_mesh;
use cortenforge::mesh::printability::{PrinterConfig, validate_for_printing};
use cortenforge::mesh::repair::validate_mesh;
use cortenforge::mesh::types::{IndexedMesh, Point3};
use cortenforge::sim::mjcf::load_model;
use cortenforge::sim::ml_chassis::{
    ActionSpace, Algorithm, EpochMetrics, LinearPolicy, ObservationSpace, TrainingBudget, VecEnv,
};
use cortenforge::sim::opt::{RicherSa, RicherSaHyperparams};
use nalgebra::{Matrix3, Vector3};
use rayon::prelude::*;

use crate::check::{GLYCEROL_MU, GLYCEROL_RHO};
use crate::helix::{Helix, V3};
use crate::rft::{advance, swim_locked};
use crate::swimmer::Swimmer;

const G: f64 = 9.81;
/// (T) the design field: enough to climb, at modest power (`coils`).
const B_DESIGN: f64 = 2e-3;
/// (m) SLA preset's least wall (PrinterConfig::sla_default), and the magnet pocket's clearance.
const SLA_WALL: f64 = 0.4e-3;
const POCKET_CLEAR: f64 = 0.05e-3;
/// (m) least gap between neighbouring turns of the wire (so they print apart and resin drains).
const GAP_MIN: f64 = 0.5e-3;
/// (m) longest swimmer: room to manoeuvre in the 50 mm jar.
const LEN_MAX: f64 = 30e-3;
/// Search box: helix radius (m), pitch angle (deg), wire radius (m).
const R_BOX: (f64, f64) = (1.0e-3, 3.0e-3);
const ANGLE_BOX: (f64, f64) = (35.0, 70.0);
const WIRE_BOX: (f64, f64) = (0.25e-3, 0.5e-3);
const TURNS: [f64; 3] = [2.0, 3.0, 4.0];
const BUDGET: usize = 60;
/// Points round the wire for the search model.
const FAST: (usize, usize) = (6, 8);

#[derive(Clone, Copy, Debug)]
struct Design {
    radius: f64,
    angle: f64,
    wire: f64,
    turns: f64,
}

impl Design {
    fn decode(a: [f64; 3], turns: f64) -> Self {
        let t = |x: f64| 0.5 * (x.clamp(-1.0, 1.0) + 1.0);
        Design {
            radius: R_BOX.0 + t(a[0]) * (R_BOX.1 - R_BOX.0),
            angle: ANGLE_BOX.0 + t(a[1]) * (ANGLE_BOX.1 - ANGLE_BOX.0),
            wire: WIRE_BOX.0 + t(a[2]) * (WIRE_BOX.1 - WIRE_BOX.0),
            turns,
        }
    }
    fn swimmer(&self) -> Swimmer {
        printable(Helix {
            radius: self.radius,
            pitch: 2.0 * PI * self.radius / self.angle.to_radians().tan(),
            turns: self.turns,
            wire: self.wire,
            hand: 1.0,
        })
    }
}

/// A helix on the printable head: the magnet plus clearance plus the SLA wall, flat open nose.
fn printable(h: Helix) -> Swimmer {
    let d = Swimmer::desk();
    Swimmer {
        helix: h,
        head_radius: d.mag_radius + POCKET_CLEAR + SLA_WALL,
        head_len: d.mag_len + POCKET_CLEAR,
        flat_nose: true,
        ..d
    }
}

#[derive(Clone, Copy, Debug)]
struct Eval {
    /// (mm/s) climb at 2 mT: top speed less sinking upright; very negative if infeasible.
    climb: f64,
    top2: f64,
    sink_up: f64,
    sink_level: f64,
    /// (Hz per mT)
    fso: f64,
    advance: f64,
    length: f64,
    gap: f64,
    mass: f64,
    feasible: bool,
}

fn evaluate(s: &Swimmer, n: (usize, usize)) -> Eval {
    let h = &s.helix;
    let gap = h.pitch - 2.0 * h.wire;
    let length = s.length();
    let short = (GAP_MIN - gap).max(0.0) + (length - LEN_MAX).max(0.0);
    let mass = s.mass();
    if short > 0.0 {
        return Eval {
            climb: -100.0 - 1e3 * short,
            top2: 0.0,
            sink_up: 0.0,
            sink_level: 0.0,
            fso: 0.0,
            advance: 0.0,
            length,
            gap,
            mass: mass.mass,
            feasible: false,
        };
    }
    let model = s.model_at(GLYCEROL_MU, n.0, n.1);
    let held = swim_locked(&model.r);
    let m = s.moment();
    let top2 = held.u.z * m * B_DESIGN / held.torque;
    let w = (mass.mass - GLYCEROL_RHO * mass.volume) * G;
    let a: Matrix3<f64> = model.r.fixed_view::<3, 3>(0, 0).into_owned();
    let a_inv = a.try_inverse().expect("drag invertible");
    let sink_up = -(a_inv * Vector3::new(0.0, 0.0, -w)).z;
    let sink_level = -(a_inv * Vector3::new(-w, 0.0, 0.0)).x;
    Eval {
        climb: (top2 - sink_up) * 1e3,
        top2,
        sink_up,
        sink_level,
        fso: m / (2.0 * PI * held.torque) * 1e-3,
        advance: advance(h, &held),
        length,
        gap,
        mass: mass.mass,
        feasible: true,
    }
}

const MJCF: &str = r#"<mujoco model="shape">
  <option gravity="0 0 0"/>
  <worldbody>
    <body name="b">
      <joint name="j" type="slide" axis="0 0 1"/>
      <geom type="sphere" size="0.01" mass="1"/>
    </body>
  </worldbody>
  <actuator>
    <motor name="radius" joint="j" gear="0" ctrllimited="true" ctrlrange="-1 1"/>
    <motor name="angle" joint="j" gear="0" ctrllimited="true" ctrlrange="-1 1"/>
    <motor name="wire" joint="j" gear="0" ctrllimited="true" ctrlrange="-1 1"/>
  </actuator>
</mujoco>"#;

type Log = Arc<Mutex<Vec<(Design, Eval)>>>;

fn make_env(log: Log, turns: f64) -> VecEnv {
    let model = Arc::new(load_model(MJCF).expect("dummy model loads"));
    let obs = ObservationSpace::builder()
        .time()
        .build(&model)
        .expect("obs space");
    let act = ActionSpace::builder()
        .ctrl(0..3)
        .build(&model)
        .expect("action space");
    VecEnv::builder(model, 1)
        .observation_space(obs)
        .action_space(act)
        .reward(move |_m, d| {
            let design = Design::decode([d.ctrl[0], d.ctrl[1], d.ctrl[2]], turns);
            let e = evaluate(&design.swimmer(), FAST);
            log.lock().expect("log lock").push((design, e));
            e.climb
        })
        .done(|_, _| true)
        .truncated(|_, _| false)
        .build()
        .expect("env builds")
}

fn describe(d: &Design) -> String {
    let h = d.swimmer().helix;
    format!(
        "R {:.2} mm, {:.1} deg (pitch {:.2} mm), wire {:.2} mm, {} turns",
        d.radius * 1e3,
        d.angle,
        h.pitch * 1e3,
        2.0 * d.wire * 1e3,
        d.turns
    )
}

fn show(label: &str, e: &Eval) {
    println!(
        "   {label:<26} climb {:+.2} mm/s = top {:.2} - sinking {:.2}; level sinking {:.2}; step-out {:.3} Hz/mT; advance {:.4}; {:.1} mm, {:.1} mg, turn gap {:.2} mm",
        e.climb,
        e.top2 * 1e3,
        e.sink_up * 1e3,
        e.sink_level * 1e3,
        e.fso,
        e.advance,
        e.length * 1e3,
        e.mass * 1e6,
        e.gap * 1e3
    );
}

/// Pull vertices onto the solid's zero set (Newton steps along the field's gradient), as in
/// stellarator `mesh`: the octree mesher's vertices sit up to most of a cell off.
fn snap(s: &Solid, m: &mut IndexedMesh, iters: usize) {
    let h = 1e-4;
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

/// The printable solid in mm, axis on z, head up: the helix swept with `pipe_spline`, a strut
/// from its end to the head, the bullet head, and the magnet pocket open at the nose.
fn solid(s: &Swimmer) -> Solid {
    let mm = 1e3;
    let h = &s.helix;
    let p = |v: V3| Point3::new(v.x * mm, v.y * mm, v.z * mm);
    let n = (32.0 * h.turns) as usize;
    let pts: Vec<Point3<f64>> = (0..=n)
        .map(|i| p(h.point(h.t_max() * i as f64 / n as f64)))
        .collect();
    let wire = Solid::pipe_spline(pts, h.wire * mm);
    let end = h.point(h.t_max());
    let base = V3::new(0.0, 0.0, s.head_z0() + s.head_radius);
    let strut = Solid::pipe_spline(vec![p(end), p((end + base) * 0.5), p(base)], h.wire * mm);
    let (r, l) = (s.head_radius * mm, s.head_len * mm);
    let z_base = base.z * mm;
    // Capsule (radius r, straight part l) centred so its straight part starts at z_base, then
    // everything above the straight part cut off: a flat nose.
    let head = Solid::capsule(r, 0.5 * l)
        .translate(Vector3::new(0.0, 0.0, z_base + 0.5 * l))
        .subtract(
            Solid::cylinder(2.0 * r, r).translate(Vector3::new(0.0, 0.0, z_base + l + r)),
        );
    let pocket_r = (s.mag_radius + POCKET_CLEAR) * mm;
    let depth = (s.mag_len + POCKET_CLEAR) * mm;
    // Pocket from z_base up through the nose (1 mm beyond it).
    let pocket = Solid::cylinder(pocket_r, 0.5 * (depth + 1.0))
        .translate(Vector3::new(0.0, 0.0, z_base + 0.5 * (depth + 1.0)));
    wire.union(strut).union(head).subtract(pocket)
}

pub fn run() {
    let t0 = Instant::now();
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("out")
        .join("shape");
    std::fs::create_dir_all(&dir).expect("create out/shape");
    let base = printable(Swimmer::desk().helix);
    println!(
        "shape: the helix for the fastest climb at {} mT, printable in resin",
        B_DESIGN * 1e3
    );
    println!(
        "   head: {:.2} mm across (magnet {:.0} mm + {:.2} mm clearance + {:.1} mm SLA wall), {:.2} mm straight, flat open nose",
        2.0 * base.head_radius * 1e3,
        2.0 * base.mag_radius * 1e3,
        POCKET_CLEAR * 1e3,
        SLA_WALL * 1e3,
        base.head_len * 1e3
    );
    println!(
        "   limits: turns at least {:.1} mm apart, at most {:.0} mm long; box R {:.1}-{:.1} mm, {:.0}-{:.0} deg, wire {:.1}-{:.1} mm",
        GAP_MIN * 1e3,
        LEN_MAX * 1e3,
        R_BOX.0 * 1e3,
        R_BOX.1 * 1e3,
        ANGLE_BOX.0,
        ANGLE_BOX.1,
        2.0 * WIRE_BOX.0 * 1e3,
        2.0 * WIRE_BOX.1 * 1e3
    );

    println!("1. The baseline (swim's helix on the printable head), search model against the full model");
    let tb = Instant::now();
    let fast = evaluate(&base, FAST);
    let tf = tb.elapsed().as_secs_f64();
    let full = evaluate(&base, (12, 16));
    show(&format!("search model ({:.1} s)", tf), &fast);
    show("full model", &full);

    println!(
        "2. Search: RicherSA, {BUDGET} evaluations for each number of turns (one seed)"
    );
    let mut csv = String::from(
        "search,turns,radius_mm,angle_deg,wire_mm,feasible,climb_mm_s,top2_mm_s,sink_up_mm_s,sink_level_mm_s,fso_hz_per_mt,advance,length_mm,mass_mg,gap_mm\n",
    );
    let mut winners = Vec::new();
    for &turns in &TURNS {
        let ts = Instant::now();
        let log: Log = Arc::new(Mutex::new(Vec::new()));
        let mut env = make_env(log.clone(), turns);
        let hp = RicherSaHyperparams {
            initial_temperature: 0.5,
            initial_proposal_std: 0.3,
            cooling_decay: 0.97,
            temperature_min: 1e-3,
            max_episode_steps: 1,
            adaptation_window: 5,
            sigma_growth_factor: 1.22,
            target_uphill_rate: 0.2,
            proposal_std_min: 1e-3,
            proposal_std_max: 1.0,
        };
        let quiet = |_: &EpochMetrics| {};
        let mut alg = RicherSa::new(Box::new(LinearPolicy::new(1, 3, &[0.0])), hp);
        alg.train(&mut env, TrainingBudget::Epochs(BUDGET), 1, &quiet);
        let rows = log.lock().expect("log lock").clone();
        let sa_secs = ts.elapsed().as_secs_f64();
        // Reference: a coarse grid over the same box, 4 x 4 x 3.
        let tg = Instant::now();
        let mut grid = Vec::new();
        for a0 in [-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0] {
            for a1 in [-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0] {
                for a2 in [-1.0, 0.0, 1.0] {
                    let d = Design::decode([a0, a1, a2], turns);
                    grid.push((d, evaluate(&d.swimmer(), FAST)));
                }
            }
        }
        let grid_secs = tg.elapsed().as_secs_f64();
        let tagged = rows
            .iter()
            .map(|r| ("richersa", r))
            .chain(grid.iter().map(|r| ("grid", r)));
        for (src, (d, e)) in tagged {
            let _ = writeln!(
                csv,
                "{src},{},{:.4},{:.3},{:.4},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.5},{:.3},{:.2},{:.3}",
                d.turns,
                d.radius * 1e3,
                d.angle,
                2.0 * d.wire * 1e3,
                e.feasible,
                e.climb,
                e.top2 * 1e3,
                e.sink_up * 1e3,
                e.sink_level * 1e3,
                e.fso,
                e.advance,
                e.length * 1e3,
                e.mass * 1e6,
                e.gap * 1e3
            );
        }
        let feasible = rows.iter().filter(|r| r.1.feasible).count();
        let top = |v: &[(Design, Eval)]| {
            *v.iter()
                .max_by(|a, b| a.1.climb.total_cmp(&b.1.climb))
                .expect("evaluations")
        };
        let (d, e) = top(&rows);
        let (gd, ge) = top(&grid);
        println!(
            "   {turns} turns, RicherSA: {} evaluations ({feasible} feasible), {sa_secs:.0} s; best {}",
            rows.len(),
            describe(&d)
        );
        show("   search model", &e);
        println!(
            "   {turns} turns, grid: {} evaluations, {grid_secs:.0} s; best {}",
            grid.len(),
            describe(&gd)
        );
        show("   search model", &ge);
        winners.push(if ge.climb > e.climb { (gd, ge) } else { (d, e) });
    }
    std::fs::write(dir.join("evals.csv"), csv).expect("write evals.csv");

    println!("3. The better of the two for each number of turns, on the full model (12 and 16 round)");
    let mut best: Option<(Design, Eval)> = None;
    for (d, _) in &winners {
        let e = evaluate(&d.swimmer(), (12, 16));
        show(&format!("{} turns", d.turns), &e);
        if best.is_none_or(|b| e.climb > b.1.climb) {
            best = Some((*d, e));
        }
    }
    let (bd, be) = best.expect("a winner");
    println!("   chosen: {}", describe(&bd));
    println!(
        "   against the baseline: climb {:+.2} -> {:+.2} mm/s; top speed at 2 mT {:.2} -> {:.2} mm/s; at 5 mT {:.2} -> {:.2} mm/s",
        full.climb,
        be.climb,
        full.top2 * 1e3,
        be.top2 * 1e3,
        full.top2 * 2.5e3,
        be.top2 * 2.5e3
    );
    println!(
        "   least field to climb: {:.2} -> {:.2} mT",
        full.sink_up / full.top2 * B_DESIGN * 1e3,
        be.sink_up / be.top2 * B_DESIGN * 1e3
    );

    println!("4. The printable part (cf-design, mm): mesh, checks, STL");
    let _ = std::fs::write(
        dir.join("best.txt"),
        format!(
            "# radius mm, pitch angle deg, wire diameter mm, turns\n{}\n{}\n{}\n{}\n",
            bd.radius * 1e3,
            bd.angle,
            2.0 * bd.wire * 1e3,
            bd.turns
        ),
    );
    println!(
        "   wrote out/shape/evals.csv and best.txt ({:.0} s in all); `shape part` makes the STL (about 15 min)",
        t0.elapsed().as_secs_f64()
    );
}

/// The printable part for a swimmer: mesh, checks, STL in `dir`. Uniform marching cubes, welded
/// at 10 nm (the octree mesher folds triangles on the curved sweep: `shape probe`), or
/// `octree` (with `snaps` Newton steps onto the surface).
fn part(sw: &Swimmer, dir: &std::path::Path, cell: Option<f64>, octree: bool, snaps: usize) {
    let t0 = Instant::now();
    let s = solid(sw);
    let cell = cell.unwrap_or((sw.helix.wire * 1e3 / 5.0).min(0.08));
    let mesh = if octree {
        let mut m = s.mesh_adaptive_par(cell).geometry;
        snap(&s, &mut m, snaps);
        m
    } else {
        let mut m = s.mesh(cell).geometry;
        cortenforge::mesh::repair::weld_vertices(&mut m, 1e-5);
        m.faces
            .retain(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
        cortenforge::mesh::repair::remove_unreferenced_vertices(&mut m);
        m
    };
    let rep = validate_mesh(&mesh);
    let vol = mesh.signed_volume();
    let model_vol = sw.mass().volume * 1e9
        - PI * (sw.mag_radius + POCKET_CLEAR).powi(2) * (sw.mag_len + POCKET_CLEAR) * 1e9;
    println!(
        "   {} at {cell:.3} mm cells: {} faces, watertight {}, manifold {}, {:.1} s",
        if octree {
            format!("octree mesh, {snaps} snap steps")
        } else {
            "marching cubes, welded".to_string()
        },
        mesh.faces.len(),
        rep.is_watertight,
        rep.is_manifold,
        t0.elapsed().as_secs_f64()
    );
    println!(
        "   volume {vol:.2} mm^3 against the model's {model_vol:.2} mm^3 less the pocket (the strut adds the rest)"
    );
    let head_z = sw.head_z0() * 1e3;
    match validate_for_printing(&mesh, &PrinterConfig::sla_default()) {
        Ok(v) => {
            let mut kinds: Vec<(String, usize, usize)> = Vec::new();
            for i in &v.issues {
                let k = format!("{:?}", i.issue_type);
                match kinds.iter_mut().find(|x| x.0 == k) {
                    Some(x) => {
                        x.1 += 1;
                        x.2 += usize::from(i.is_critical());
                    }
                    None => kinds.push((k, 1, usize::from(i.is_critical()))),
                }
            }
            println!(
                "   SLA check, printed upright (head up): {}",
                if v.is_printable() {
                    "printable"
                } else {
                    "not printable as is"
                }
            );
            for (k, n, c) in &kinds {
                println!("      {k}: {n} ({c} critical)");
            }
            let mut tw: Vec<_> = v.thin_walls.iter().collect();
            tw.sort_by(|a, b| a.thickness.total_cmp(&b.thickness));
            let in_head = tw.iter().filter(|t| t.center.z >= head_z).count();
            println!(
                "      thin-wall regions: {} ({} in the head, z >= {head_z:.2} mm); thinnest:",
                tw.len(),
                in_head
            );
            for t in tw.iter().take(5) {
                println!(
                    "         {:.3} mm at ({:.2}, {:.2}, {:.2}), {} faces, {:.3} mm^2",
                    t.thickness,
                    t.center.x,
                    t.center.y,
                    t.center.z,
                    t.faces.len(),
                    t.area
                );
            }
            for i in v
                .issues
                .iter()
                .filter(|i| format!("{:?}", i.issue_type) == "SelfIntersecting")
            {
                println!(
                    "      self-intersection: {} at {:?}",
                    i.description,
                    i.location.map(|l| (l.x, l.y, l.z))
                );
            }
            if let Some(i) = v
                .issues
                .iter()
                .filter(|i| i.is_critical() && format!("{:?}", i.issue_type) == "ExcessiveOverhang")
                .max_by(|a, b| a.affected_elements.len().cmp(&b.affected_elements.len()))
            {
                println!("      largest critical overhang: {}", i.description);
            }
        }
        Err(e) => println!("   printability check failed: {e}"),
    }
    let path = dir.join("swimmer.stl");
    save_mesh(&mesh, &path).expect("write STL");
    println!("   wrote {}", path.display());
}

/// `shape part [cell mm] [octree [snap steps]]`: the part alone, for the design in
/// out/shape/best.txt; marching cubes unless `octree` is given.
pub fn run_part() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("out")
        .join("shape");
    let txt = std::fs::read_to_string(dir.join("best.txt")).expect("out/shape/best.txt (run `shape` first)");
    let v: Vec<f64> = txt
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter_map(|l| l.split_whitespace().next()?.parse().ok())
        .collect();
    let d = Design {
        radius: v[0] * 1e-3,
        angle: v[1],
        wire: 0.5 * v[2] * 1e-3,
        turns: v[3],
    };
    println!("shape part: {}", describe(&d));
    let args: Vec<String> = std::env::args().skip(3).collect();
    let cell = args.first().and_then(|a| a.parse().ok());
    let octree = args.get(1).is_some_and(|a| a == "octree");
    let snaps = args.get(2).and_then(|a| a.parse().ok()).unwrap_or(2);
    part(&d.swimmer(), &dir, cell, octree, snaps);
}

/// One mesh's defects under the SLA check: (overlapping pairs reported, critical thin walls,
/// thinnest wall mm, watertight and manifold).
fn defects(m: &IndexedMesh) -> (String, usize, f64, bool) {
    let rep = validate_mesh(m);
    match validate_for_printing(m, &PrinterConfig::sla_default()) {
        Ok(v) => {
            let si = v
                .issues
                .iter()
                .find(|i| format!("{:?}", i.issue_type) == "SelfIntersecting")
                .map(|i| i.description.split_whitespace().next().unwrap_or("?").to_string())
                .unwrap_or_else(|| "0".into());
            let crit_thin = v.thin_walls.iter().filter(|t| t.thickness < 0.2).count();
            let thinnest = v
                .thin_walls
                .iter()
                .map(|t| t.thickness)
                .fold(f64::INFINITY, f64::min);
            (si, crit_thin, thinnest, rep.is_watertight && rep.is_manifold)
        }
        Err(_) => ("check failed".into(), 0, f64::NAN, false),
    }
}

/// `shape probe`: where the part's mesh defects come from. Three pieces of the chosen design,
/// each meshed by the octree mesher (raw and snapped) and by uniform marching cubes (welded).
pub fn run_probe() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("out")
        .join("shape");
    let txt = std::fs::read_to_string(dir.join("best.txt")).expect("out/shape/best.txt");
    let v: Vec<f64> = txt
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter_map(|l| l.split_whitespace().next()?.parse().ok())
        .collect();
    let d = Design {
        radius: v[0] * 1e-3,
        angle: v[1],
        wire: 0.5 * v[2] * 1e-3,
        turns: v[3],
    };
    let sw = d.swimmer();
    let mm = 1e3;
    let h = &sw.helix;
    let b = h.wire * mm;
    let cell = (b / 5.0).min(0.08);
    let n = (32.0 * h.turns) as usize;
    let pts: Vec<Point3<f64>> = (0..=n)
        .map(|i| {
            let q = h.point(h.t_max() * i as f64 / n as f64);
            Point3::new(q.x * mm, q.y * mm, q.z * mm)
        })
        .collect();
    let straight: Vec<Point3<f64>> = (0..=n)
        .map(|i| Point3::new(10.0 * i as f64 / n as f64, 0.0, 0.0))
        .collect();
    let pieces: Vec<(&str, Solid)> = vec![
        ("capsule, 10 mm, wire radius", Solid::capsule(b, 5.0)),
        ("straight pipe_spline, 10 mm", Solid::pipe_spline(straight, b)),
        ("the helix wire (pipe_spline)", Solid::pipe_spline(pts, b)),
        ("the whole part", solid(&sw)),
    ];
    println!(
        "shape probe: {} at {cell:.3} mm cells; SLA check (thin = under 0.2 mm, critical)",
        describe(&d)
    );
    println!(
        "   {:<30} {:<22} {:>8} {:>10} {:>12} {:>9} {:>7}",
        "piece", "mesher", "faces", "overlaps", "thin (crit)", "thinnest", "secs"
    );
    for (name, s) in &pieces {
        for mesher in ["octree", "octree + snap 2", "marching cubes + weld"] {
            let t0 = Instant::now();
            let mut m = match mesher {
                "marching cubes + weld" => {
                    let mut m = s.mesh(cell).geometry;
                    cortenforge::mesh::repair::weld_vertices(&mut m, 1e-5);
                    m.faces
                        .retain(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
                    cortenforge::mesh::repair::remove_unreferenced_vertices(&mut m);
                    m
                }
                _ => s.mesh_adaptive_par(cell).geometry,
            };
            if mesher == "octree + snap 2" {
                snap(s, &mut m, 2);
            }
            let secs = t0.elapsed().as_secs_f64();
            let (si, thin, thinnest, ok) = defects(&m);
            println!(
                "   {name:<30} {mesher:<22} {:>8} {si:>10} {thin:>12} {thinnest:>9.3} {secs:>7.1}{}",
                m.faces.len(),
                if ok { "" } else { "  (not watertight/manifold)" }
            );
        }
    }
}
