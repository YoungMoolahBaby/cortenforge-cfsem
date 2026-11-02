//! `layout`: widening the pad's sweet spot.
//!
//! The receiver is Qi receiver example 1 (free space, as in `pad`). A pad is one or more round
//! single-layer coils of the A11 Litz; with several coils the pad drives the one that couples best
//! (as Qi's free-positioning designs do), the others left open. Score: the area of receiver
//! positions (centre within |x| <= 40 mm, |y| <= 24 mm, long side along x) where the link's best
//! efficiency (both sides tuned, best load, 100 kHz) is at least 85%.
//!
//! Every round coil is a sum of rings, so the receiver's coupling to one ring is tabulated once
//! over radius and offset, for the two layer heights (exact circles: the receiver's vector
//! potential from cfsem, integrated round the ring by the trapezoid rule), and interpolated
//! (Catmull-Rom in all three). Coil inductance: cfsem's coaxial ring mutuals plus Wien's ring
//! formula for each turn. Benchmarks: the A11 coil, and Qi's free-positioning design A6 (three
//! overlapping rectangular coils; Part 4 v1.2.3, section 2.2.6) computed directly.

use std::f64::consts::PI;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_circular_filament;
use cfsem::physics::linear_filament::vector_potential_linear_filament;
use cortenforge::sim::mjcf::load_model;
use cortenforge::sim::ml_chassis::{
    ActionSpace, Algorithm, LinearPolicy, ObservationSpace, TrainingBudget, VecEnv,
};
use cortenforge::sim::opt::{RicherSa, RicherSaHyperparams, Sa, SaHyperparams};
use cortenforge::sim::rl::{Cem, CemHyperparams};
use rayon::prelude::*;

use crate::loops::{Loop, mutual};
use crate::pad::{F, LITZ_AC, RX_WIRE, RX_Z, TX_WIRE, eta_max, l_coil, rx_coil, tx_coil};

const W: f64 = 2.0 * PI * F;
/// Litz of the A11 / A6 coils: 105 x 0.08 mm, 1.15 mm across (the least turn pitch).
const PITCH: f64 = 1.15e-3;
fn ohm_per_m() -> f64 {
    1.72e-8 / (105.0 * PI / 4.0 * 0.08e-3 * 0.08e-3) * LITZ_AC
}
/// Layers: the top layer's centre RX_Z under the receiver's; the next layer one pitch lower.
const SEPS: [f64; 2] = [RX_Z, RX_Z + PITCH];
const ETA: f64 = 0.85;
/// Receiver positions: x -40..40 mm, y 0..24 mm (mirrored in y), 2 mm steps.
const XS: (i32, i32) = (-20, 20);
const YS: i32 = 12;
const STEP: f64 = 2e-3;
/// Pad footprint for arrays: turn centrelines within |x| <= 52 mm.
const HALF_PAD: f64 = 0.052;
/// Least pitch of the array: the outer pair (same layer) keeps its outermost turns a wire apart.
fn s_min(r_out: f64) -> f64 {
    r_out + 0.5 * PITCH
}

/// The receiver as straight segments (each side one exact segment), 1 A per turn.
struct Segs {
    s: (Vec<f64>, Vec<f64>, Vec<f64>),
    d: (Vec<f64>, Vec<f64>, Vec<f64>),
    i: Vec<f64>,
    r: Vec<f64>,
}

fn rx_segs() -> Segs {
    let (mut s, mut d) = ((vec![], vec![], vec![]), (vec![], vec![], vec![]));
    for t in rx_coil(10.0) {
        let (a, b) = t.segments();
        for (p, q) in a.iter().zip(&b) {
            s.0.push(p.0);
            s.1.push(p.1);
            s.2.push(p.2);
            d.0.push(q.0);
            d.1.push(q.1);
            d.2.push(q.2);
        }
    }
    let n = s.0.len();
    Segs {
        s,
        d,
        i: vec![1.0; n],
        r: vec![RX_WIRE; n],
    }
}

/// (H) the receiver's coupling to one exact ring (radius r, centre (cx, cy), `sep` below the
/// receiver's plane): A of the receiver round the ring, trapezoid rule at about `h` spacing.
fn ring_m(rx: &Segs, r: f64, cx: f64, cy: f64, sep: f64, h: f64) -> f64 {
    let n = ((2.0 * PI * r / h).ceil() as usize).max(48);
    let th: Vec<f64> = (0..n).map(|k| 2.0 * PI * k as f64 / n as f64).collect();
    let px: Vec<f64> = th.iter().map(|t| cx + r * t.cos()).collect();
    let py: Vec<f64> = th.iter().map(|t| cy + r * t.sin()).collect();
    let pz = vec![-sep; n];
    let (mut ax, mut ay, mut az) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    vector_potential_linear_filament(
        (&px, &py, &pz),
        (&rx.s.0, &rx.s.1, &rx.s.2),
        (&rx.d.0, &rx.d.1, &rx.d.2),
        &rx.i,
        &rx.r,
        (&mut ax, &mut ay, &mut az),
    )
    .expect("lengths");
    th.iter()
        .enumerate()
        .map(|(k, t)| -ax[k] * t.sin() + ay[k] * t.cos())
        .sum::<f64>()
        * 2.0
        * PI
        * r
        / n as f64
}

/// M(ring radius, |dx|, |dy|) per layer, on a grid.
struct Table {
    r0: f64,
    dr: f64,
    nr: usize,
    nx: usize,
    ny: usize,
    data: Vec<f64>,
}

impl Table {
    fn build(rx: &Segs) -> Self {
        let (r0, dr, nr, nx, ny) = (3e-3, 1e-3, 38, 46, 21); // r 3..40 mm; dx 0..90, dy 0..40 mm
        let idx: Vec<(usize, usize, usize, usize)> = (0..2)
            .flat_map(|l| {
                (0..nr).flat_map(move |ir| {
                    (0..nx).flat_map(move |ix| (0..ny).map(move |iy| (l, ir, ix, iy)))
                })
            })
            .collect();
        let data = idx
            .par_iter()
            .map(|&(l, ir, ix, iy)| {
                ring_m(
                    rx,
                    r0 + dr * ir as f64,
                    STEP * ix as f64,
                    STEP * iy as f64,
                    SEPS[l],
                    1e-3,
                )
            })
            .collect();
        Table {
            r0,
            dr,
            nr,
            nx,
            ny,
            data,
        }
    }

    fn at(&self, l: usize, ir: i64, ix: i64, iy: i64) -> f64 {
        let ir = ir.clamp(0, self.nr as i64 - 1) as usize;
        let ix = (ix.unsigned_abs() as usize).min(self.nx - 1);
        let iy = (iy.unsigned_abs() as usize).min(self.ny - 1);
        self.data[((l * self.nr + ir) * self.nx + ix) * self.ny + iy]
    }

    /// Catmull-Rom in radius and both offsets (the offsets mirror at zero).
    fn m(&self, l: usize, r: f64, dx: f64, dy: f64) -> f64 {
        let w = |t: f64| {
            [
                0.5 * (-t * t * t + 2.0 * t * t - t),
                0.5 * (3.0 * t * t * t - 5.0 * t * t + 2.0),
                0.5 * (-3.0 * t * t * t + 4.0 * t * t + t),
                0.5 * (t * t * t - t * t),
            ]
        };
        let split = |u: f64| (u.floor() as i64, u - u.floor());
        let (ir, tr) = split((r - self.r0) / self.dr);
        let (ix, tx) = split(dx.abs() / STEP);
        let (iy, ty) = split(dy.abs() / STEP);
        let (wr, wx, wy) = (w(tr), w(tx), w(ty));
        let off = |i: i64, k: usize| i + k as i64 - 1;
        let mut s = 0.0;
        for (a, wa) in wr.iter().enumerate() {
            for (b, wb) in wx.iter().enumerate() {
                for (c, wc) in wy.iter().enumerate() {
                    s += wa * wb * wc * self.at(l, off(ir, a), off(ix, b), off(iy, c));
                }
            }
        }
        s
    }
}

/// (H) a single-layer round coil's inductance: Wien's formula per turn (wire radius 0.575 mm)
/// plus cfsem's coaxial ring mutuals.
fn l_rings(rings: &[f64]) -> f64 {
    let a = TX_WIRE;
    let wien = |r: f64| {
        MU_0 * r
            * ((1.0 + a * a / (8.0 * r * r)) * (8.0 * r / a).ln() - 0.0083 * a * a / (r * r) - 1.75)
    };
    let n = rings.len();
    let mut l: f64 = rings.iter().map(|&r| wien(r)).sum();
    let mut out = vec![0.0; n];
    for i in 0..n {
        flux_circular_filament(
            (&[rings[i]], &[0.0], &[1.0]),
            &[a],
            (rings, &vec![0.0; n]),
            &mut out,
        )
        .expect("lengths");
        l += out
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, m)| m)
            .sum::<f64>();
    }
    l
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Design {
    r_in: f64,
    r_out: f64,
    n: usize,
    /// Array pitch (0 for a single coil): coils at -s (lower layer), 0 (top), +s (lower layer).
    s: f64,
}

impl Design {
    fn rings(&self) -> Vec<f64> {
        (0..self.n)
            .map(|k| self.r_in + (self.r_out - self.r_in) * k as f64 / (self.n - 1) as f64)
            .collect()
    }
    fn coils(&self) -> Vec<(f64, usize)> {
        if self.s == 0.0 {
            vec![(0.0, 0)]
        } else {
            vec![(-self.s, 1), (0.0, 0), (self.s, 1)]
        }
    }
    fn describe(&self) -> String {
        format!(
            "{} coil{} r {:.1}-{:.1} mm, {} turns (pitch {:.2} mm){}",
            self.coils().len(),
            if self.s == 0.0 { "" } else { "s" },
            self.r_in * 1e3,
            self.r_out * 1e3,
            self.n,
            (self.r_out - self.r_in) / (self.n - 1) as f64 * 1e3,
            if self.s == 0.0 {
                String::new()
            } else {
                format!(", {:.1} mm apart", self.s * 1e3)
            }
        )
    }
    /// From actions in [-1, 1]; three: an array (r_out <= 26 mm so the outer pair fits apart).
    fn decode(a: &[f64], three: bool) -> Option<Self> {
        let t = |i: usize| 0.5 * (a[i].clamp(-1.0, 1.0) + 1.0);
        let r_in = 3e-3 + t(0) * 27e-3;
        let r_max = if three { 26e-3 } else { 40e-3 };
        if r_max - r_in < 2.0 * PITCH {
            return None;
        }
        let r_out = r_in + 2.0 * PITCH + t(1) * (r_max - r_in - 2.0 * PITCH);
        let n_max = ((r_out - r_in) / PITCH).floor() as usize + 1;
        let n = 2 + ((t(2) * (n_max - 1) as f64) as usize).min(n_max - 2);
        let s = if three {
            let lo = s_min(r_out);
            lo + t(3) * (HALF_PAD - r_out - lo).max(0.0)
        } else {
            0.0
        };
        if three && s + r_out > HALF_PAD + 1e-9 {
            return None;
        }
        Some(Design { r_in, r_out, n, s })
    }
}

struct Scorer {
    table: Table,
    l2: f64,
    q2: f64,
}

/// Receiver positions and their area weights (y = 0 counts half, then the y mirror doubles).
fn positions() -> Vec<(f64, f64, f64)> {
    (XS.0..=XS.1)
        .flat_map(|i| {
            (0..=YS).map(move |j| {
                (
                    STEP * i as f64,
                    STEP * j as f64,
                    if j == 0 { 1.0 } else { 2.0 } * STEP * STEP * 1e6,
                )
            })
        })
        .collect()
}

impl Scorer {
    /// (soft area, hard area in mm^2, centred k, Q1) for coupling sums from `m_at`.
    fn area(&self, l1: f64, r1: f64, m_at: &dyn Fn(f64, f64) -> f64) -> (f64, f64, f64, f64) {
        let q1 = W * l1 / r1;
        let (mut soft, mut hard) = (0.0, 0.0);
        for (x, y, w) in positions() {
            let k = m_at(x, y).abs() / (l1 * self.l2).sqrt();
            let e = eta_max(k, q1, self.q2);
            soft += w / (1.0 + (-(e - ETA) / 0.004).exp());
            if e >= ETA {
                hard += w;
            }
        }
        (soft, hard, m_at(0.0, 0.0).abs() / (l1 * self.l2).sqrt(), q1)
    }

    fn score(&self, d: &Design) -> (f64, f64, f64, f64) {
        let rings = d.rings();
        let l1 = l_rings(&rings);
        let r1 = rings.iter().map(|r| 2.0 * PI * r).sum::<f64>() * ohm_per_m();
        let coils = d.coils();
        let m_at = |x: f64, y: f64| {
            coils
                .iter()
                .map(|&(cx, l)| {
                    rings
                        .iter()
                        .map(|&r| self.table.m(l, r, x - cx, y))
                        .sum::<f64>()
                        .abs()
                })
                .fold(0.0, f64::max)
        };
        self.area(l1, r1, &m_at)
    }
}

/// A pad scored without the table (exact rings at every position): (area with efficiency >= 85%
/// in mm^2, soft area, coil inductance, centred k, Q1). `s` = 0 for one coil.
pub(crate) fn direct_area(r_in: f64, r_out: f64, n: usize, s: f64) -> (f64, f64, f64, f64, f64) {
    let d = Design { r_in, r_out, n, s };
    let rx = rx_segs();
    let rx_turns = rx_coil(1e-3);
    let l2 = l_coil(&rx_turns, RX_WIRE, None, None);
    let r2 = rx_turns.iter().map(Loop::length).sum::<f64>() * 1.72e-8
        / (26.0 * PI / 4.0 * 0.08e-3 * 0.08e-3)
        * LITZ_AC;
    let (rings, coils) = (d.rings(), d.coils());
    let l1 = l_rings(&rings);
    let q1 = W * l1 / (rings.iter().map(|r| 2.0 * PI * r).sum::<f64>() * ohm_per_m());
    let k_at = |x: f64, y: f64| {
        coils
            .iter()
            .map(|&(cx, l)| {
                rings
                    .iter()
                    .map(|&r| ring_m(&rx, r, cx - x, -y, SEPS[l], 1e-3))
                    .sum::<f64>()
                    .abs()
            })
            .fold(0.0, f64::max)
            / (l1 * l2).sqrt()
    };
    let (hard, soft) = positions()
        .par_iter()
        .map(|&(x, y, w)| {
            let e = eta_max(k_at(x, y), q1, W * l2 / r2);
            (
                if e >= ETA { w } else { 0.0 },
                w / (1.0 + (-(e - ETA) / 0.004).exp()),
            )
        })
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
    (hard, soft, l1, k_at(0.0, 0.0), q1)
}

/// One-step episode for sim-opt / sim-rl: actions carry the design, the reward is its soft area.
fn mjcf(n: usize) -> String {
    let motors: String = (0..n)
        .map(|i| {
            format!(
                r#"<motor name="a{i}" joint="j" gear="0" ctrllimited="true" ctrlrange="-1 1"/>"#
            )
        })
        .collect();
    format!(
        r#"<mujoco model="pad"><option gravity="0 0 0"/><worldbody><body name="b"><joint name="j" type="slide" axis="0 0 1"/>
<geom type="sphere" size="0.01" mass="1"/></body></worldbody><actuator>{motors}</actuator></mujoco>"#
    )
}

type Log = Arc<Mutex<Vec<f64>>>;

fn make_env(sc: Arc<Scorer>, n_envs: usize, n_act: usize, three: bool, log: Log) -> VecEnv {
    let model = Arc::new(load_model(&mjcf(n_act)).expect("dummy model loads"));
    let obs = ObservationSpace::builder()
        .time()
        .build(&model)
        .expect("obs space");
    let act = ActionSpace::builder()
        .ctrl(0..n_act)
        .build(&model)
        .expect("action space");
    VecEnv::builder(model, n_envs)
        .observation_space(obs)
        .action_space(act)
        .reward(move |_m, d| {
            let a: Vec<f64> = (0..n_act).map(|i| d.ctrl[i]).collect();
            let s = Design::decode(&a, three)
                .map(|d| sc.score(&d).0)
                .unwrap_or(0.0);
            log.lock().expect("log").push(s);
            s
        })
        .done(|_, _| true)
        .truncated(|_, _| false)
        .build()
        .expect("env builds")
}

/// Runs one optimizer; returns (shipped design's hard area, best soft area seen, evaluations, s).
fn optimize(
    sc: &Arc<Scorer>,
    name: &str,
    three: bool,
    budget: usize,
    seed: u64,
) -> (Option<Design>, f64, usize, f64) {
    let n_act = if three { 4 } else { 3 };
    let log: Log = Arc::new(Mutex::new(Vec::new()));
    let t0 = Instant::now();
    let quiet = |_: &cortenforge::sim::ml_chassis::EpochMetrics| {};
    let policy = || Box::new(LinearPolicy::new(1, n_act, &[0.0]));
    let params = match name {
        "SA" => {
            let mut env = make_env(sc.clone(), 1, n_act, three, log.clone());
            let hp = SaHyperparams {
                initial_temperature: 50.0,
                proposal_std: 0.1,
                cooling_decay: 0.99,
                temperature_min: 0.1,
                max_episode_steps: 1,
            };
            let mut a = Sa::new(policy(), hp);
            a.train(&mut env, TrainingBudget::Epochs(budget), seed, &quiet);
            a.best_artifact().params
        }
        "RicherSA" => {
            let mut env = make_env(sc.clone(), 1, n_act, three, log.clone());
            let hp = RicherSaHyperparams {
                initial_temperature: 50.0,
                initial_proposal_std: 0.3,
                cooling_decay: 0.99,
                temperature_min: 0.1,
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
        _ => {
            let pop = 16;
            let mut env = make_env(sc.clone(), pop, n_act, three, log.clone());
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
    };
    let secs = t0.elapsed().as_secs_f64();
    let log = log.lock().expect("log");
    let best = log.iter().cloned().fold(0.0, f64::max);
    let a: Vec<f64> = (0..n_act)
        .map(|i| f64::from(params[n_act + i].tanh() as f32))
        .collect();
    (Design::decode(&a, three), best, log.len(), secs)
}

pub fn run() {
    let t0 = Instant::now();
    let rx = rx_segs();
    let rx_turns = rx_coil(1e-3);
    let l2 = l_coil(&rx_turns, RX_WIRE, None, None);
    let r2 = rx_turns.iter().map(Loop::length).sum::<f64>() * 1.72e-8
        / (26.0 * PI / 4.0 * 0.08e-3 * 0.08e-3)
        * LITZ_AC;
    println!(
        "Pad layout: receiver example 1 in free space (L2 {:.2} uH, Q2 {:.0}); score = area of receiver centres in \
         |x| <= 40, |y| <= 24 mm with best efficiency >= {:.0}%",
        l2 * 1e6,
        W * l2 / r2,
        ETA * 100.0
    );

    // 1. Checks of the ring model.
    let a11 = tx_coil(256);
    let a11_rings: Vec<f64> = (0..10)
        .map(|k| 10.825e-3 + (21.425e-3 - 10.825e-3) * k as f64 / 9.0)
        .collect();
    println!(
        "1. A11 inductance: rings (Wien + cfsem coaxial mutuals) {:.4} uH, polygons (256-gons) {:.4} uH",
        l_rings(&a11_rings) * 1e6,
        l_coil(&a11, TX_WIRE, None, None) * 1e6
    );
    let one = ring_m(&rx, 15.3e-3, 7.1e-3, 3.3e-3, RX_Z, 1e-3);
    let half = ring_m(&rx, 15.3e-3, 7.1e-3, 3.3e-3, RX_Z, 0.5e-3);
    println!(
        "   one ring (r 15.3 mm, offset 7.1, 3.3 mm): trapezoid at 1 mm {one:.6e} H, at 0.5 mm {half:.6e} H ({:+.1e})",
        one / half - 1.0
    );
    let t1 = Instant::now();
    let table = Table::build(&rx);
    println!(
        "   ring table: {} entries in {:.1} s",
        table.data.len(),
        t1.elapsed().as_secs_f64()
    );
    for (r, x, y, l) in [
        (15.3e-3, 7.1e-3, 3.3e-3, 0),
        (8.7e-3, 31.9e-3, 17.3e-3, 1),
        (33.4e-3, 0.9e-3, 22.1e-3, 0),
    ] {
        let want = ring_m(&rx, r, -x, y, SEPS[l], 1e-3);
        let got = table.m(l, r, x, y);
        println!(
            "   interpolated vs direct, r {:.1} mm at ({:.1}, {:.1}) mm, layer {l}: {got:.6e} vs {want:.6e} H ({:+.1e} of the centred value)",
            r * 1e3,
            x * 1e3,
            y * 1e3,
            (got - want) / ring_m(&rx, r, 0.0, 0.0, SEPS[l], 1e-3)
        );
    }
    let sc = Arc::new(Scorer {
        table,
        l2,
        q2: W * l2 / r2,
    });
    // The A11 via the table against `pad`'s polygon map at its centre and (12, 8) mm.
    for (x, y) in [(0.0, 0.0), (12e-3, 8e-3)] {
        let got: f64 = a11_rings.iter().map(|&r| sc.table.m(0, r, x, y)).sum();
        let rxp: Vec<Loop> = rx_coil(2e-3)
            .iter()
            .map(|t| t.placed(0.0, (x, y, RX_Z)))
            .collect();
        let want: f64 = tx_coil(256)
            .iter()
            .map(|t| rxp.iter().map(|r| mutual(t, r)).sum::<f64>())
            .sum();
        println!(
            "   A11 coupling at ({:.0}, {:.0}) mm: table {got:.6e}, 256-gons {want:.6e} H ({:+.1e})",
            x * 1e3,
            y * 1e3,
            got / want - 1.0
        );
    }

    // 2. Benchmarks.
    let a11d = Design {
        r_in: 10.825e-3,
        r_out: 21.425e-3,
        n: 10,
        s: 0.0,
    };
    let s11 = sc.score(&a11d);
    println!("2. Benchmarks:");
    println!(
        "   Qi A11 ({}): {:.0} mm^2 (centred k {:.3}, Q1 {:.0})",
        a11d.describe(),
        s11.1,
        s11.2,
        s11.3
    );
    // A6: rectangular coils 45.2 x 53.2 mm outside (the array along their 45.2 mm width), 19.5 x
    // 27.5 inside, 12 turns, 1.5 mm thick; outer pair 49.2 mm apart in the lower layer, the middle
    // one on top. Its receiver gap: 0.75 + 2.0 + 2.5 + 0.3 mm.
    let p6 = (45.2e-3 - 19.5e-3) / 2.0 / 12.0;
    let a6_turns = |piece: f64| -> Vec<Loop> {
        (0..12)
            .map(|k| {
                Loop::rect_by(
                    0.0,
                    0.0,
                    0.0,
                    45.2e-3 - (2 * k + 1) as f64 * p6,
                    53.2e-3 - (2 * k + 1) as f64 * p6,
                    piece,
                    1.0,
                )
            })
            .collect()
    };
    let a6 = a6_turns(2e-3);
    let l6 = l_coil(&a6_turns(1e-3), TX_WIRE, None, None);
    let l6_top = l_coil(
        &a6_turns(1e-3),
        TX_WIRE,
        Some(-1.5e-3 - 0.75e-3 - 0.15e-3),
        None,
    );
    let l6_bot = l_coil(&a6_turns(1e-3), TX_WIRE, Some(-0.75e-3 - 0.15e-3), None);
    println!(
        "   A6 coil inductance: free {:.2} uH; over an ideal ferrite plane {:.2} uH (top coil), {:.2} uH (lower coil); spec 11.5 / 12.5 uH +- 10%",
        l6 * 1e6,
        l6_top * 1e6,
        l6_bot * 1e6
    );
    let r6 = a6.iter().map(Loop::length).sum::<f64>() * ohm_per_m();
    let z6 = 0.75e-3 + 2.0e-3 + 2.5e-3 + 0.3e-3;
    let rx_src = rx_coil(10.0);
    let pos = positions();
    let m6: Vec<f64> = pos
        .par_iter()
        .map(|&(x, y, _)| {
            let r: Vec<Loop> = rx_src.iter().map(|t| t.placed(0.0, (x, y, 0.0))).collect();
            [
                (-24.6e-3, -z6 - 1.5e-3),
                (0.0, -z6),
                (24.6e-3, -z6 - 1.5e-3),
            ]
            .iter()
            .map(|&(cx, cz)| {
                let tx: Vec<Loop> = a6.iter().map(|t| t.placed(0.0, (cx, 0.0, cz))).collect();
                r.iter()
                    .map(|s| tx.iter().map(|t| mutual(s, t)).sum::<f64>())
                    .sum::<f64>()
                    .abs()
            })
            .fold(0.0, f64::max)
        })
        .collect();
    let lookup = |x: f64, y: f64| {
        m6[pos
            .iter()
            .position(|p| (p.0 - x).abs() < 1e-9 && (p.1 - y).abs() < 1e-9)
            .expect("grid")]
    };
    let s6 = sc.area(l6, r6, &lookup);
    println!(
        "   Qi A6 (three 45.2 x 53.2 mm coils, 24.6 mm apart, 12 turns): {:.0} mm^2 (centred k {:.3}, Q1 {:.0})",
        s6.1, s6.2, s6.3
    );

    // 3. Grid searches.
    let t2 = Instant::now();
    let mut cands: Vec<Design> = Vec::new();
    for i in 0..=27 {
        let r_in = 3e-3 + 1e-3 * i as f64;
        for o in 0..=37 {
            let r_out = r_in + 2.0 * PITCH + 1e-3 * o as f64;
            if r_out > 40e-3 + 1e-12 {
                break;
            }
            let n_max = ((r_out - r_in) / PITCH).floor() as usize + 1;
            for n in 2..=n_max {
                cands.push(Design {
                    r_in,
                    r_out,
                    n,
                    s: 0.0,
                });
            }
        }
    }
    let single: Vec<(Design, (f64, f64, f64, f64))> =
        cands.par_iter().map(|d| (*d, sc.score(d))).collect();
    let best1 = single
        .iter()
        .cloned()
        .fold((cands[0], (0.0, 0.0, 0.0, 0.0)), |a, b| {
            if b.1.0 > a.1.0 { b } else { a }
        });
    println!(
        "3. Grid, one coil (r_in 3-30 mm, r_out to 40 mm in 1 mm steps, every turn count): {} designs in {:.1} s; best {}: {:.0} mm^2 \
         (centred k {:.3}, Q1 {:.0})",
        single.len(),
        t2.elapsed().as_secs_f64(),
        best1.0.describe(),
        best1.1.1,
        best1.1.2,
        best1.1.3
    );
    let t3 = Instant::now();
    let mut c3: Vec<Design> = Vec::new();
    for i in 0..=22 {
        let r_in = 3e-3 + 1e-3 * i as f64;
        for o in 0..=23 {
            let r_out = r_in + 2.0 * PITCH + 1e-3 * o as f64;
            if r_out > 26e-3 + 1e-12 {
                break;
            }
            let n_max = ((r_out - r_in) / PITCH).floor() as usize + 1;
            for n in 2..=n_max {
                let mut s = s_min(r_out);
                while s + r_out <= HALF_PAD + 1e-9 {
                    c3.push(Design { r_in, r_out, n, s });
                    s += 2e-3;
                }
            }
        }
    }
    let three: Vec<(Design, (f64, f64, f64, f64))> =
        c3.par_iter().map(|d| (*d, sc.score(d))).collect();
    let best3 = three
        .iter()
        .cloned()
        .fold((c3[0], (0.0, 0.0, 0.0, 0.0)), |a, b| {
            if b.1.0 > a.1.0 { b } else { a }
        });
    println!(
        "   grid, three overlapping coils (r_out to 26 mm, pitch 2 mm steps): {} designs in {:.1} s; best {}: {:.0} mm^2 \
         (centred k {:.3}, Q1 {:.0})",
        three.len(),
        t3.elapsed().as_secs_f64(),
        best3.0.describe(),
        best3.1.1,
        best3.1.2,
        best3.1.3
    );

    // Refine each grid best on a fine local grid: radii +-1.5 mm in 0.1 mm steps, turns +-2, the
    // array pitch +-2 mm in 0.2 mm steps.
    let refine = |d: Design| -> (Design, (f64, f64, f64, f64)) {
        let mut v = Vec::new();
        for i in -15..=15 {
            for o in -15..=15 {
                let (r_in, r_out) = (d.r_in + 1e-4 * i as f64, d.r_out + 1e-4 * o as f64);
                if r_in < 3e-3 || r_out > (if d.s == 0.0 { 40e-3 } else { 26e-3 }) {
                    continue;
                }
                for dn in -2i64..=2 {
                    let n = (d.n as i64 + dn).max(2) as usize;
                    if (r_out - r_in) / ((n - 1) as f64) < PITCH - 1e-12 {
                        continue;
                    }
                    let ss: Vec<f64> = if d.s == 0.0 {
                        vec![0.0]
                    } else {
                        (-10..=10).map(|k| d.s + 2e-4 * k as f64).collect()
                    };
                    for s in ss {
                        if s != 0.0 && (s < s_min(r_out) - 1e-12 || s + r_out > HALF_PAD + 1e-9) {
                            continue;
                        }
                        v.push(Design { r_in, r_out, n, s });
                    }
                }
            }
        }
        v.par_iter().map(|d| (*d, sc.score(d))).reduce(
            || (d, sc.score(&d)),
            |a, b| if b.1.0 > a.1.0 { b } else { a },
        )
    };
    let t4 = Instant::now();
    let (fine1, fine3) = (refine(best1.0), refine(best3.0));
    println!(
        "   refined (0.1 mm radii, 0.2 mm pitch, turns +-2; {:.1} s): one coil {}: {:.0} mm^2 (soft {:.1}); three coils {}: {:.0} mm^2 (soft {:.1})",
        t4.elapsed().as_secs_f64(),
        fine1.0.describe(),
        fine1.1.1,
        fine1.1.0,
        fine3.0.describe(),
        fine3.1.1,
        fine3.1.0
    );

    // 4. Direct check of the winners (no table: exact rings at every position).
    for (name, d) in [("one coil", fine1.0), ("three coils", fine3.0)] {
        let rings = d.rings();
        let coils = d.coils();
        let ms: Vec<f64> = pos
            .par_iter()
            .map(|&(x, y, _)| {
                coils
                    .iter()
                    .map(|&(cx, l)| {
                        rings
                            .iter()
                            .map(|&r| ring_m(&rx, r, cx - x, -y, SEPS[l], 1e-3))
                            .sum::<f64>()
                            .abs()
                    })
                    .fold(0.0, f64::max)
            })
            .collect();
        let look = |x: f64, y: f64| {
            ms[pos
                .iter()
                .position(|p| (p.0 - x).abs() < 1e-9 && (p.1 - y).abs() < 1e-9)
                .expect("grid")]
        };
        let l1 = l_rings(&rings);
        let r1 = rings.iter().map(|r| 2.0 * PI * r).sum::<f64>() * ohm_per_m();
        let direct = sc.area(l1, r1, &look);
        println!(
            "4. {name}, direct (no table): {:.0} mm^2, soft {:.1} vs table {:.1} mm^2; inductance {:.2} uH",
            direct.1,
            direct.0,
            sc.score(&d).0,
            l1 * 1e6
        );
    }

    // 5. CortenForge's optimizers against the grid optima (soft area).
    println!(
        "5. sim-opt (SA, RicherSA) and sim-rl (CEM) as one-step episodes, 3 seeds: the soft area of the design each ships \
         (its best artifact) and the best it saw, against the refined grid"
    );
    for (label, three, best) in [("one coil", false, fine1.1), ("three coils", true, fine3.1)] {
        let mut top: Option<(Design, f64)> = None;
        for budget in [400, 1600] {
            for name in ["SA", "RicherSA", "CEM"] {
                let runs: Vec<(Option<Design>, f64, usize, f64)> = (0..3)
                    .map(|seed| optimize(&sc, name, three, budget, seed))
                    .collect();
                let shipped: Vec<String> = runs
                    .iter()
                    .map(|r| match r.0 {
                        Some(d) => {
                            let s = sc.score(&d).0;
                            if top.is_none_or(|t| s > t.1) {
                                top = Some((d, s));
                            }
                            format!("{s:.1}")
                        }
                        None => "invalid".into(),
                    })
                    .collect();
                let seen: Vec<String> = runs.iter().map(|r| format!("{:.1}", r.1)).collect();
                println!(
                    "   {label:<11} {name:<8} budget {budget:>4}: shipped {}, best seen {} (refined grid {:.1}); {:.1} s a run",
                    shipped.join(" / "),
                    seen.join(" / "),
                    best.0,
                    runs.iter().map(|r| r.3).sum::<f64>() / 3.0
                );
            }
        }
        if let Some((d, s)) = top {
            let full = sc.score(&d);
            println!(
                "   {label}: best shipped {}: soft {s:.1}, {:.0} mm^2; L {:.2} uH, centred k {:.3}, Q1 {:.0}",
                d.describe(),
                full.1,
                l_rings(&d.rings()) * 1e6,
                full.2,
                full.3
            );
        }
    }
    println!("   total {:.0} s", t0.elapsed().as_secs_f64());
}
