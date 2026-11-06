//! `optimize`: the coil search. sim-rl's CEM moves the 45 coil coefficients (`search::Design`)
//! to raise `search::evaluate`'s score.
//!
//! CortenForge's optimizers only run on a `VecEnv` built from a sim-core model (see
//! docs/CORTENFORGE.md, sim-opt), so as in coil-designer the design rides on a dummy MJCF model:
//! one inert body, one zero-gear motor per coefficient, a `LinearPolicy` whose biases carry the
//! design (its observation weight is scaled to zero), a one-step episode whose reward reads
//! `ctrl` and scores the coils. CEM's `best_artifact` returns an elite mean it never evaluated
//! (docs/CORTENFORGE.md, sim-rl), so every evaluation is logged here and the best evaluated
//! design is what gets saved.
//!
//! Writes `out/opt/evals.csv` (one row per evaluation) and `out/opt/best.txt` (the best design's
//! coefficients, rewritten whenever the best improves).

use std::fmt::Write as _;
use std::io::Write as _;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use cortenforge::sim::mjcf::load_model;
use cortenforge::sim::ml_chassis::{
    ActionSpace, Algorithm, EpochMetrics, LinearPolicy, ObservationSpace, TrainingBudget, VecEnv,
};
use cortenforge::sim::rl::{Cem, CemHyperparams};

use crate::design;
use crate::search::{self, Design, N_PARAMS, Score};

fn mjcf() -> String {
    let mut s = String::from(
        r#"<mujoco model="coils">
  <option gravity="0 0 0"/>
  <worldbody>
    <body name="b">
      <joint name="j" type="slide" axis="0 0 1"/>
      <geom type="sphere" size="0.01" mass="1"/>
    </body>
  </worldbody>
  <actuator>
"#,
    );
    for k in 0..N_PARAMS {
        let _ = writeln!(
            s,
            r#"    <motor name="c{k}" joint="j" gear="0" ctrllimited="true" ctrlrange="-1 1"/>"#
        );
    }
    s.push_str("  </actuator>\n</mujoco>\n");
    s
}

struct Log {
    rows: Vec<(Design, Score)>,
    best: Option<(Design, Score)>,
    csv: std::fs::File,
    dir: std::path::PathBuf,
    t0: Instant,
}

impl Log {
    fn push(&mut self, d: Design, s: Score) {
        let n = self.rows.len() + 1;
        let mut line = format!(
            "{n},{:.4},{:.1},{:.5},{:.5},{:.2},{:.2},{:.2},{},{:.2}",
            s.value, s.r_good, s.iota_axis, s.iota_min, s.rc_min, s.gap, s.penalty, s.note, s.secs
        );
        for x in &d.x {
            let _ = write!(line, ",{x:.6}");
        }
        let _ = writeln!(self.csv, "{line}");
        let better = self.best.as_ref().is_none_or(|b| s.value > b.1.value);
        if better {
            let mut txt = format!(
                "# value {:.4}, r_good {:.1} mm, iota axis {:.5}, iota min {:.5}, rc_min {:.2} mm, gap {:.2} mm, evaluation {n}\n",
                s.value, s.r_good, s.iota_axis, s.iota_min, s.rc_min, s.gap
            );
            for x in &d.x {
                let _ = writeln!(txt, "{x:.9}");
            }
            std::fs::write(self.dir.join("best.txt"), txt).expect("write best.txt");
            println!(
                "  [{n:>4}, {:>5.0} s] new best {:.2}: r_good {:.0} mm, iota axis {:.4}, min {:.4}, rc {:.1} mm, gap {:.1} mm",
                self.t0.elapsed().as_secs_f64(),
                s.value,
                s.r_good,
                s.iota_axis,
                s.iota_min,
                s.rc_min,
                s.gap
            );
            self.best = Some((d.clone(), s.clone()));
        }
        self.rows.push((d, s));
    }
}

fn make_env(pop: usize, log: Arc<Mutex<Log>>) -> VecEnv {
    let model = Arc::new(load_model(&mjcf()).expect("dummy model loads"));
    let obs = ObservationSpace::builder().time().build(&model).expect("obs space");
    let act = ActionSpace::builder()
        .ctrl(0..N_PARAMS)
        .build(&model)
        .expect("action space");
    let base = design::baseline();
    VecEnv::builder(model, pop)
        .observation_space(obs)
        .action_space(act)
        .reward(move |_m, d| {
            let x: Vec<f64> = (0..N_PARAMS).map(|k| d.ctrl[k]).collect();
            let design = Design { x };
            let s = search::evaluate(&design, &base);
            let v = s.value;
            log.lock().expect("log lock").push(design, s);
            v
        })
        .done(|_, _| true)
        .truncated(|_, _| false)
        .build()
        .expect("env builds")
}

/// `optimize [generations] [population] [seed]`.
pub fn run(args: &[String]) {
    let gens: usize = args.first().and_then(|a| a.parse().ok()).unwrap_or(24);
    let pop: usize = args.get(1).and_then(|a| a.parse().ok()).unwrap_or(24);
    let seed: u64 = args.get(2).and_then(|a| a.parse().ok()).unwrap_or(1);
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out").join("opt");
    std::fs::create_dir_all(&dir).expect("create out/opt");
    let mut csv = std::fs::File::create(dir.join("evals.csv")).expect("evals.csv");
    let mut head = String::from("n,value,r_good_mm,iota_axis,iota_min,rc_min_mm,gap_mm,penalty,note,secs");
    for k in 0..N_PARAMS {
        let _ = write!(head, ",x{k}");
    }
    let _ = writeln!(csv, "{head}");
    println!(
        "coil search: CEM, {gens} generations x {pop} candidates (seed {seed}), {N_PARAMS} coefficients; \
         score = r_good + {} x min(iota, {}) - penalties",
        search::IOTA_WEIGHT,
        search::IOTA_CAP
    );
    let log = Arc::new(Mutex::new(Log {
        rows: Vec::new(),
        best: None,
        csv,
        dir: dir.clone(),
        t0: Instant::now(),
    }));

    // The baseline first, as the reference every candidate is measured against.
    let base = design::baseline();
    let b = search::evaluate(&Design::zero(), &base);
    println!(
        "baseline: {:.2} (r_good {:.0} mm, iota axis {:.4}, min {:.4}), {:.1} s",
        b.value, b.r_good, b.iota_axis, b.iota_min, b.secs
    );
    log.lock().expect("log lock").push(Design::zero(), b);

    let mut env = make_env(pop, log.clone());
    let hp = CemHyperparams {
        elite_fraction: 0.25,
        noise_std: 0.2,
        noise_decay: 0.93,
        noise_min: 0.02,
        max_episode_steps: 1,
    };
    // One observation (time) scaled to zero: params are N_PARAMS dead weights then the biases,
    // which start at zero (the baseline).
    let mut cem = Cem::new(Box::new(LinearPolicy::new(1, N_PARAMS, &[0.0])), hp);
    let progress = {
        let log = log.clone();
        move |m: &EpochMetrics| {
            let l = log.lock().expect("log lock");
            let recent = &l.rows[l.rows.len().saturating_sub(pop)..];
            let mean = recent.iter().map(|r| r.1.value).sum::<f64>() / recent.len().max(1) as f64;
            let best = l.best.as_ref().map_or(f64::NAN, |b| b.1.value);
            println!(
                "generation {:>3}: mean {mean:.2}, best so far {best:.2} ({} evaluations, {:.0} s)",
                m.epoch + 1,
                l.rows.len(),
                l.t0.elapsed().as_secs_f64()
            );
        }
    };
    cem.train(&mut env, TrainingBudget::Epochs(gens), seed, &progress);

    let l = log.lock().expect("log lock");
    if let Some((_, s)) = &l.best {
        println!(
            "done: best {:.2} (r_good {:.0} mm, iota axis {:.4}, min {:.4}, rc {:.1} mm, gap {:.1} mm) after {} evaluations, {:.0} s; saved to {}",
            s.value,
            s.r_good,
            s.iota_axis,
            s.iota_min,
            s.rc_min,
            s.gap,
            l.rows.len(),
            l.t0.elapsed().as_secs_f64(),
            dir.join("best.txt").display()
        );
    }
    let _ = std::io::stdout().flush();
}

/// Read a saved design (`best.txt` format) back: `path`, or `out/opt/best.txt` by default.
pub fn load_best(path: Option<&str>) -> Design {
    let p = match path {
        Some(p) => std::path::PathBuf::from(p),
        None => std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out").join("opt").join("best.txt"),
    };
    let txt = std::fs::read_to_string(&p).expect("out/opt/best.txt (run `optimize` first)");
    let x: Vec<f64> = txt
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.trim().parse().expect("a coefficient"))
        .collect();
    assert_eq!(x.len(), N_PARAMS, "best.txt holds {} coefficients", x.len());
    Design { x }
}

/// `best`: the saved best design at full resolution, beside the baseline: the search's own
/// score again, then `trace` (128 filaments a coil, 600 periods a line) and `mesh`.
pub fn best(args: &[String]) {
    let base = design::baseline();
    let d = load_best(args.iter().find(|a| !a.starts_with("--")).map(String::as_str));
    let s = search::evaluate(&d, &base);
    println!(
        "best design, search score again: {:.2} (r_good {:.0} mm, iota axis {:.4}, min {:.4}, rc {:.1} mm, gap {:.1} mm)",
        s.value, s.r_good, s.iota_axis, s.iota_min, s.rc_min, s.gap
    );
    let coils = d.coils(&base);
    let free = d.free_coils(&base);
    for (j, c) in free.iter().enumerate() {
        let b = base.coil(j);
        println!(
            "  free coil {j}: length {:.1} mm (baseline {:.1}), out of plane by up to {:.1} mm",
            c.length(256) * 1e3,
            b.length(256) * 1e3,
            out_of_plane(c, base.phi(j)) * 1e3
        );
    }
    if !args.iter().any(|a| a == "--no-trace") {
        crate::poincare::trace_coils(&base, &coils, "optimized", "poincare_best.svg");
    }
    if !args.iter().any(|a| a == "--no-mesh") {
        crate::meshes::mesh_coils(&base, &coils, "mesh_best");
    }
}

/// (m) largest distance of a coil from the poloidal plane at `phi`.
fn out_of_plane(c: &crate::coils::FourierCurve, phi: f64) -> f64 {
    let n = [-phi.sin(), phi.cos(), 0.0];
    c.sample(256)
        .iter()
        .map(|p| (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]).abs())
        .fold(0.0, f64::max)
}
