//! `optimize`: Helmholtz spacing and pack shape for the largest uniform sphere.
//!
//! CortenForge's gradient-free optimizers (sim-opt `Sa`, `RicherSa`, `Pt`; sim-rl `Cem`) only run
//! on a `VecEnv`, which needs a sim-core `Model`. So the design problem is dressed as a one-step
//! episode: a dummy MJCF model with two inert actuators, a `LinearPolicy` whose biases are the
//! design parameters (its obs weight is scaled to zero), and a reward that reads `ctrl`, builds
//! the coil and returns its uniform-sphere radius. A brute-force grid is the reference.

use std::sync::{Arc, Mutex};
use std::time::Instant;

use cortenforge::sim::mjcf::load_model;
use cortenforge::sim::ml_chassis::{
    ActionSpace, Algorithm, LinearPolicy, ObservationSpace, TrainingBudget, VecEnv,
};
use cortenforge::sim::opt::{Pt, PtHyperparams, RicherSa, RicherSaHyperparams, Sa, SaHyperparams};
use cortenforge::sim::rl::{Cem, CemHyperparams};

use crate::coil::{Coil, Wire, uniform_radius};
use crate::pack::helmholtz;

const R_MEAN: f64 = 0.05; // m
const TURNS: usize = 100; // per coil, approximately (layers x round(TURNS / layers))
const EPS: f64 = 1e-4; // uniformity target
const LAYERS: (usize, usize) = (1, 25); // radial layers, inclusive
const SPACING: (f64, f64) = (0.9, 1.6); // spacing / R_MEAN
/// (m) least axial gap between the two packs: room for two 2.5 mm former walls (an assumption).
const GAP_MIN: f64 = 5e-3;
const BUDGETS: [usize; 2] = [400, 1600]; // reward evaluations per run
const SEEDS: u64 = 3;

/// A candidate: spacing as a fraction of the mean radius, and radial layer count.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Design {
    spacing: f64,
    layers: usize,
}

impl Design {
    /// From two actions in [-1, 1], in equal bins for the layer count; `fixed` pins the layers.
    fn decode(a0: f64, a1: f64, fixed: Option<usize>) -> Self {
        let t0 = 0.5 * (a0.clamp(-1.0, 1.0) + 1.0);
        let t1 = 0.5 * (a1.clamp(-1.0, 1.0) + 1.0);
        let spacing = SPACING.0 + t0 * (SPACING.1 - SPACING.0);
        let n = LAYERS.1 - LAYERS.0 + 1;
        let layers = fixed.unwrap_or(LAYERS.0 + ((t1 * n as f64) as usize).min(n - 1));
        Design { spacing, layers }
    }

    fn per_layer(&self) -> usize {
        (TURNS as f64 / self.layers as f64).round() as usize
    }

    fn coil(&self) -> Coil {
        helmholtz(
            R_MEAN,
            self.spacing * R_MEAN,
            self.layers,
            self.per_layer(),
            Wire::awg(22),
        )
    }

    /// (m) axial gap between the packs.
    fn gap(&self) -> f64 {
        self.spacing * R_MEAN - self.coil().packs[0].width()
    }

    /// (mm) radius of the EPS-uniform sphere; 0 if the packs are closer than [GAP_MIN].
    fn score(&self) -> f64 {
        if self.gap() < GAP_MIN {
            return 0.0;
        }
        let c = self.coil();
        let f = |r: &[f64], z: &[f64]| c.field(1.0, r, z);
        uniform_radius(&f, EPS, 0.3 * R_MEAN) * 1e3
    }

    fn describe(&self) -> String {
        let p = self.coil().packs[0];
        format!(
            "spacing {:.4} R, {} x {} turns ({:.1} mm deep x {:.1} mm wide)",
            self.spacing,
            self.layers,
            self.per_layer(),
            p.depth() * 1e3,
            p.width() * 1e3
        )
    }
}

/// The dummy model: one inert body and two zero-gear motors whose `ctrl` carries the design.
const MJCF: &str = r#"<mujoco model="design">
  <option gravity="0 0 0"/>
  <worldbody>
    <body name="b">
      <joint name="j" type="slide" axis="0 0 1"/>
      <geom type="sphere" size="0.01" mass="1"/>
    </body>
  </worldbody>
  <actuator>
    <motor name="spacing" joint="j" gear="0" ctrllimited="true" ctrlrange="-1 1"/>
    <motor name="shape" joint="j" gear="0" ctrllimited="true" ctrlrange="-1 1"/>
  </actuator>
</mujoco>"#;

/// Every reward evaluation: (design, score in mm).
type Log = Arc<Mutex<Vec<(Design, f64)>>>;

fn make_env(n_envs: usize, log: Log, fixed: Option<usize>) -> VecEnv {
    let model = Arc::new(load_model(MJCF).expect("dummy model loads"));
    let obs = ObservationSpace::builder()
        .time()
        .build(&model)
        .expect("obs space");
    let act = ActionSpace::builder()
        .ctrl(0..2)
        .build(&model)
        .expect("action space");
    VecEnv::builder(model, n_envs)
        .observation_space(obs)
        .action_space(act)
        .reward(move |_m, d| {
            let design = Design::decode(d.ctrl[0], d.ctrl[1], fixed);
            let s = design.score();
            log.lock().expect("log lock").push((design, s));
            s
        })
        .done(|_, _| true)
        .truncated(|_, _| false)
        .build()
        .expect("env builds")
}

/// One obs (time) scaled to zero, two actions: params are [w0, w1, b0, b1].
fn policy() -> Box<LinearPolicy> {
    Box::new(LinearPolicy::new(1, 2, &[0.0]))
}

/// The design a user would ship: the algorithm's best artifact, decoded the way the env does
/// (tanh, then f32 at the action boundary).
fn shipped(params: &[f64], fixed: Option<usize>) -> Design {
    let a = |b: f64| f64::from(b.tanh() as f32);
    Design::decode(a(params[2]), a(params[3]), fixed)
}

struct Run {
    shipped: Design,
    shipped_score: f64,
    best_seen: f64,
    evals: usize,
    secs: f64,
}

fn run_one(name: &str, seed: u64, budget: usize, fixed: Option<usize>) -> Run {
    let log: Log = Arc::new(Mutex::new(Vec::new()));
    let t0 = Instant::now();
    let quiet = |_: &cortenforge::sim::ml_chassis::EpochMetrics| {};
    let params = match name {
        "SA" => {
            let mut env = make_env(1, log.clone(), fixed);
            let hp = SaHyperparams {
                initial_temperature: 0.5,
                proposal_std: 0.1,
                cooling_decay: 0.99,
                temperature_min: 1e-3,
                max_episode_steps: 1,
            };
            let mut a = Sa::new(policy(), hp);
            a.train(&mut env, TrainingBudget::Epochs(budget), seed, &quiet);
            a.best_artifact().params
        }
        "RicherSA" => {
            let mut env = make_env(1, log.clone(), fixed);
            let hp = RicherSaHyperparams {
                initial_temperature: 0.5,
                initial_proposal_std: 0.3,
                cooling_decay: 0.99,
                temperature_min: 1e-3,
                max_episode_steps: 1,
                adaptation_window: 5,
                sigma_growth_factor: 1.22,
                target_uphill_rate: 0.2,
                proposal_std_min: 1e-4,
                proposal_std_max: 1.0,
            };
            let mut a = RicherSa::new(policy(), hp);
            a.train(&mut env, TrainingBudget::Epochs(budget), seed, &quiet);
            a.best_artifact().params
        }
        "PT" => {
            let mut env = make_env(1, log.clone(), fixed);
            let k = 4;
            let hp = PtHyperparams {
                n_chains: k,
                t_min: 0.005,
                t_max: 0.5,
                proposal_std: 0.1,
                max_episode_steps: 1,
                swap_interval_epochs: 1,
            };
            let mut a = Pt::new(policy(), hp);
            a.train(&mut env, TrainingBudget::Epochs(budget / k), seed, &quiet);
            a.best_artifact().params
        }
        "CEM" => {
            let pop = 16;
            let mut env = make_env(pop, log.clone(), fixed);
            let hp = CemHyperparams {
                elite_fraction: 0.25,
                noise_std: 0.3,
                noise_decay: 0.9,
                noise_min: 0.005,
                max_episode_steps: 1,
            };
            let mut a = Cem::new(policy(), hp);
            a.train(&mut env, TrainingBudget::Epochs(budget / pop), seed, &quiet);
            a.best_artifact().params
        }
        _ => unreachable!(),
    };
    let secs = t0.elapsed().as_secs_f64();
    let log = log.lock().expect("log lock");
    let best_seen = log.iter().map(|x| x.1).fold(0.0, f64::max);
    let d = shipped(&params, fixed);
    Run {
        shipped: d,
        shipped_score: d.score(),
        best_seen,
        evals: log.len(),
        secs,
    }
}

pub fn run() {
    println!(
        "Helmholtz pair, mean radius {:.0} mm, ~{TURNS} turns/coil of AWG 22: largest sphere with \
         |B - B0|/B0 <= {EPS:.0e}",
        R_MEAN * 1e3
    );
    println!(
        "Search space: spacing {:.2}-{:.2} R, {}-{} layers (turns per layer = round({TURNS} / layers))",
        SPACING.0, SPACING.1, LAYERS.0, LAYERS.1
    );

    // Reference: brute force, spacing step 0.002 R for every layer count, then 5e-5 R refinement.
    let t0 = Instant::now();
    let mut best = (
        Design {
            spacing: 1.0,
            layers: 10,
        },
        0.0,
    );
    let mut per_layers = Vec::new();
    let steps = ((SPACING.1 - SPACING.0) / 0.002).round() as usize;
    let mut evals = 0;
    for layers in LAYERS.0..=LAYERS.1 {
        let mut b = (0.0, 0.0);
        for k in 0..=steps {
            let d = Design {
                spacing: SPACING.0 + 0.002 * k as f64,
                layers,
            };
            let s = d.score();
            evals += 1;
            if s > b.1 {
                b = (d.spacing, s);
            }
        }
        if b.1 > 0.0 {
            for k in -40..=40 {
                let d = Design {
                    spacing: b.0 + 5e-5 * k as f64,
                    layers,
                };
                let s = d.score();
                evals += 1;
                if s > b.1 {
                    b = (d.spacing, s);
                }
            }
        }
        per_layers.push((layers, b));
        if b.1 > best.1 {
            best = (
                Design {
                    spacing: b.0,
                    layers,
                },
                b.1,
            );
        }
    }
    println!(
        "Brute force: {evals} evaluations in {:.1} s ({:.1} ms each); best spacing per layer count",
        t0.elapsed().as_secs_f64(),
        t0.elapsed().as_secs_f64() * 1e3 / evals as f64
    );
    for (layers, (sp, s)) in &per_layers {
        let d = Design {
            spacing: *sp,
            layers: *layers,
        };
        if *s > 0.0 {
            let c = d.coil();
            let b0 = c.field(1.0, &[0.0], &[0.0]).1[0];
            let p1 = (1e-3 / b0).powi(2) * c.resistance();
            println!(
                "  {:<62} gap {:5.1} mm  {s:.3} mm  {:.3} mT/A  {p1:.2} W for 1 mT",
                d.describe(),
                d.gap() * 1e3,
                b0 * 1e3
            );
        } else {
            println!(
                "  {layers} layers: no spacing in range leaves a {:.0} mm gap",
                GAP_MIN * 1e3
            );
        }
    }
    println!("  optimum: {} -> {:.3} mm", best.0.describe(), best.1);

    // The ideal thin pair with the same turns, spacing optimized the same way.
    let n = TURNS as f64;
    let thin_at = |sp: f64| {
        let loops = [
            (R_MEAN, -0.5 * sp * R_MEAN, n),
            (R_MEAN, 0.5 * sp * R_MEAN, n),
        ];
        let f = |r: &[f64], z: &[f64]| crate::pack::thin_field(&loops, r, z);
        uniform_radius(&f, EPS, 0.3 * R_MEAN) * 1e3
    };
    let mut tb = (1.0, thin_at(1.0));
    for k in 0..=2000 {
        let sp = 0.95 + 1e-4 * k as f64;
        let s = thin_at(sp);
        if s > tb.1 {
            tb = (sp, s);
        }
    }
    println!(
        "Thin loops for comparison: {:.3} mm at spacing 1.0000 R, {:.3} mm at the best spacing {:.4} R",
        thin_at(1.0),
        tb.1,
        tb.0
    );

    // Both parameters at two budgets, then spacing alone at a fixed 10 x 10 pack.
    let ten = per_layers[10 - LAYERS.0].1.1;
    let mut cases: Vec<(String, usize, Option<usize>, f64)> = BUDGETS
        .iter()
        .map(|&b| {
            (
                format!("spacing and layers, {b} evaluations"),
                b,
                None,
                best.1,
            )
        })
        .collect();
    cases.push((
        "spacing only at 10 x 10, 200 evaluations".into(),
        200,
        Some(10),
        ten,
    ));
    for (label, budget, fixed, reference) in cases {
        println!(
            "CortenForge optimizers: {label}, {SEEDS} seeds, reference {reference:.3} mm \n             (shipped = best_artifact)"
        );
        for name in ["SA", "RicherSA", "PT", "CEM"] {
            let mut gaps = Vec::new();
            for seed in 0..SEEDS {
                let r = run_one(name, seed, budget, fixed);
                let gap = r.shipped_score / reference - 1.0;
                gaps.push(gap);
                println!(
                    "  {name:<8} seed {seed}: shipped {:.3} mm ({gap:+.2}%), best seen {:.3} mm, \
                     {} evals, {:.1} s; {}",
                    r.shipped_score,
                    r.best_seen,
                    r.evals,
                    r.secs,
                    r.shipped.describe(),
                    gap = gap * 100.0
                );
            }
            gaps.sort_by(f64::total_cmp);
            println!(
                "  {name:<8} median gap {:+.2}%, worst {:+.2}%",
                gaps[gaps.len() / 2] * 100.0,
                gaps[0] * 100.0
            );
        }
    }
}
