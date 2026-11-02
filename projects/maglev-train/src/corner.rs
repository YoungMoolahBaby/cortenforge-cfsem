//! `corner`: one corner of the cart in CortenForge's simulator (sim-core), with the magnetics
//! from cfsem.
//!
//! The corner: a 12 x 4 mm N42 disc with a control coil wound round it, under the steel rail.
//! Upward pull F(g, I) = a(g) + b(g) I + c(g) I^2, tabulated from cfsem (magnet and coil against
//! their images). The coil links the magnet's flux Psi(g), so its voltage is
//! V = R I + L(g) dI/dt + (I L'(g) + Psi'(g)) dg/dt. Energy says b = -dPsi/dg and c = -L'/2,
//! both checked below.
//! The body is a sim-core slide joint; the force goes in through `qfrc_applied` each step and the
//! coil current is integrated alongside.

use std::f64::consts::PI;
use std::fmt::Write as _;
use std::sync::Arc;

use cortenforge::sim::core::{Data, Model};
use cortenforge::sim::mjcf::load_model;

use crate::em::{Coil, Wire, force_lorentz, images, magnet_loops, mutual, self_inductance};

const MASS: f64 = 0.5 / 4.0; // kg a corner
const G: f64 = 9.81;
const SUPPLY: f64 = 12.0; // V
const NOISE: f64 = 1e-6; // m, gap sensor noise (1 sigma); the sweep in step 4 shows why
const SAMPLE: f64 = 2000.0; // Hz, gap sensor and outer loop
const I_MAX: f64 = 10.0; // A, current command limit
const RAMP: f64 = 0.3; // s, lift-off setpoint ramp

/// Tabulated corner magnetics on a uniform gap grid.
struct Table {
    g0: f64,
    dg: f64,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    l: Vec<f64>,
    psi: Vec<f64>,
    r: f64,
}

impl Table {
    fn build(coil: &Coil) -> Self {
        let (g0, g1, n) = (0.002, 0.012, 201);
        let dg = (g1 - g0) / (n - 1) as f64;
        let l0 = self_inductance(coil);
        let (mut a, mut b, mut c, mut l, mut psi) = (vec![], vec![], vec![], vec![], vec![]);
        for k in 0..n {
            let g = g0 + k as f64 * dg;
            let (m, im) = magnet_loops(0.006, 0.004, 1.3, 24, g);
            let t = coil.turns_at(g);
            let (mi, ti) = (images(&m), images(&t));
            a.push(force_lorentz(&m, im, &mi, im));
            b.push(force_lorentz(&m, im, &ti, 1.0) + force_lorentz(&t, 1.0, &mi, im));
            c.push(force_lorentz(&t, 1.0, &ti, 1.0));
            l.push(l0 + mutual(&t, &ti));
            let both: Vec<(f64, f64)> = m.iter().copied().chain(mi).collect();
            psi.push(mutual(&t, &both) * im);
        }
        Table {
            g0,
            dg,
            a,
            b,
            c,
            l,
            psi,
            r: coil.resistance(),
        }
    }
    /// Catmull-Rom value and derivative of a column at gap g.
    fn eval(&self, col: &[f64], g: f64) -> (f64, f64) {
        let x = ((g - self.g0) / self.dg).clamp(1.0, (col.len() - 3) as f64);
        let i = x.floor() as usize;
        let t = x - i as f64;
        let (p0, p1, p2, p3) = (col[i - 1], col[i], col[i + 1], col[i + 2]);
        let v = 0.5
            * (2.0 * p1
                + (-p0 + p2) * t
                + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t * t
                + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t * t * t);
        let d = 0.5
            * ((-p0 + p2)
                + 2.0 * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t
                + 3.0 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t * t)
            / self.dg;
        (v, d)
    }
    fn force(&self, g: f64, i: f64) -> f64 {
        self.eval(&self.a, g).0 + self.eval(&self.b, g).0 * i + self.eval(&self.c, g).0 * i * i
    }
    /// Gap where the magnet alone carries the corner.
    fn zero_power_gap(&self, weight: f64) -> f64 {
        let (mut lo, mut hi) = (
            self.g0 + 2.0 * self.dg,
            self.g0 + (self.a.len() - 3) as f64 * self.dg,
        );
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.eval(&self.a, mid).0 > weight {
                lo = mid
            } else {
                hi = mid
            }
        }
        0.5 * (lo + hi)
    }
}

fn corner_model(dt: f64) -> Arc<Model> {
    let xml = format!(
        r#"<mujoco model="corner">
  <option timestep="{dt}" gravity="0 0 -{G}"/>
  <worldbody>
    <body name="corner" pos="0 0 0">
      <joint name="z" type="slide" axis="0 0 1"/>
      <geom type="box" size="0.01 0.01 0.005" mass="{MASS}"/>
    </body>
  </worldbody>
</mujoco>"#
    );
    Arc::new(load_model(&xml).expect("corner model loads"))
}

/// The controllers: cascade of a gap PD (sampled, noisy sensor) giving a current command, an
/// analog-like PI current loop giving the coil voltage (clamped to the supply), and optionally a
/// slow zero-power loop moving the gap setpoint until the mean current is zero.
#[derive(Clone, Copy)]
struct Gains {
    /// (m) gap sensor noise, 1 sigma.
    noise: f64,
    /// (Hz) corner of the first-order filter on the measured gap rate.
    d_filter: f64,
    kp: f64,
    kd: f64,
    kpi: f64,
    kii: f64,
    zero_power: f64,
}

#[derive(Clone, Copy, PartialEq)]
enum Drive {
    Open,
    Shorted,
    Control(bool),
}

struct Log {
    t: Vec<f64>,
    gap: Vec<f64>,
    current: Vec<f64>,
    power: Vec<f64>,
}

/// Run one corner for `t_end` s from gap `start`. `payload` adds mass at `t_payload`.
#[allow(clippy::too_many_arguments)]
fn simulate(
    tab: &Table,
    drive: Drive,
    gains: Gains,
    start: f64,
    g_ref0: f64,
    t_end: f64,
    dt: f64,
    payload: (f64, f64),
) -> Log {
    let model = corner_model(dt);
    let mut data: Data = model.make_data();
    // qpos is the upward displacement from the start; gap = start - qpos.
    let mut i = 0.0;
    let mut i_int = 0.0;
    let mut g_ref = g_ref0;
    let mut i_cmd: f64 = 0.0;
    let mut last_meas = start;
    let mut dmeas_f = 0.0;
    let steps_per_sample = ((1.0 / SAMPLE) / dt).round() as usize;
    let mut rng = 0x2545_f491_4f6c_dd1d_u64;
    let mut gauss = || {
        // Box-Muller from a xorshift.
        let mut u = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            (rng >> 11) as f64 / (1u64 << 53) as f64
        };
        let (a, b) = (u().max(1e-300), u());
        (-2.0 * a.ln()).sqrt() * (2.0 * PI * b).cos()
    };
    let mut log = Log {
        t: vec![],
        gap: vec![],
        current: vec![],
        power: vec![],
    };
    let n = (t_end / dt).round() as usize;
    for step in 0..n {
        let t = step as f64 * dt;
        let gap = start - data.qpos[0];
        let gdot = -data.qvel[0];
        if !(0.0015..0.02).contains(&gap) {
            break; // stuck to the rail or dropped off the table
        }
        let extra = if t >= payload.0 { payload.1 } else { 0.0 };
        // Controller.
        let v = match drive {
            Drive::Open | Drive::Shorted => 0.0,
            Drive::Control(zp) => {
                if step % steps_per_sample == 0 {
                    let meas = gap + gains.noise * gauss();
                    let ts = 1.0 / SAMPLE;
                    // Derivative through a first-order filter.
                    let raw = (meas - last_meas) / ts;
                    let alpha = 1.0 - (-2.0 * PI * gains.d_filter * ts).exp();
                    dmeas_f += alpha * (raw - dmeas_f);
                    last_meas = meas;
                    if zp && t > RAMP {
                        // Mean current above zero: the magnet is short of the load, so move
                        // closer (smaller gap), where it pulls harder.
                        g_ref -= gains.zero_power * i_cmd * ts;
                    }
                    // Lift-off: ramp the setpoint from the starting gap over RAMP seconds.
                    let target = if t < RAMP {
                        start + (g_ref - start) * t / RAMP
                    } else {
                        g_ref
                    };
                    i_cmd = (gains.kp * (meas - target) + gains.kd * dmeas_f).clamp(-I_MAX, I_MAX);
                }
                let e = i_cmd - i;
                let v_raw = gains.kpi * e + gains.kii * i_int;
                // Anti-windup: integrate only while the voltage is not pinned at the supply.
                if v_raw.abs() < SUPPLY || v_raw.signum() != e.signum() {
                    i_int += e * dt;
                }
                v_raw.clamp(-SUPPLY, SUPPLY)
            }
        };
        // Coil: V = R I + L dI/dt + (I L' + Psi') dg/dt (open coil: I stays 0).
        if drive != Drive::Open {
            let (l, dl) = tab.eval(&tab.l, gap);
            let (_, dpsi) = tab.eval(&tab.psi, gap);
            let di = (v - tab.r * i - (i * dl + dpsi) * gdot) / l;
            i += di * dt;
        }
        let f = tab.force(gap, i);
        data.qfrc_applied[0] = f - extra * G;
        if step % 10 == 0 {
            log.t.push(t);
            log.gap.push(gap);
            log.current.push(i);
            log.power.push(i * i * tab.r);
        }
        data.step(&model).expect("step");
    }
    log
}

/// Growth rate from rest: a released offset grows as d0 cosh(p t), so p = acosh(10) / t10 with
/// t10 the first time |gap - g0| reaches ten times the start (interpolated between log samples).
fn growth_rate(log: &Log, g0: f64) -> f64 {
    let d0 = (log.gap[0] - g0).abs();
    for k in 1..log.t.len() {
        let (a, b) = ((log.gap[k - 1] - g0).abs(), (log.gap[k] - g0).abs());
        if b >= 10.0 * d0 {
            let t = log.t[k - 1] + (10.0 * d0 - a) / (b - a) * (log.t[k] - log.t[k - 1]);
            return 10f64.acosh() / t;
        }
    }
    f64::NAN
}

/// The largest real root of L s^3 + R s^2 + (b^2/m - k L/m) s - k R / m = 0 (shorted coil).
fn shorted_root(l: f64, r: f64, b: f64, k: f64, m: f64) -> f64 {
    let f = |s: f64| l * s * s * s + r * s * s + (b * b / m - k * l / m) * s - k * r / m;
    let (mut lo, mut hi) = (0.0, (k / m).sqrt() * 2.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 { lo = mid } else { hi = mid }
    }
    0.5 * (lo + hi)
}

pub fn run() {
    let coil = Coil {
        r_in: 0.0065,
        layers: 8,
        per_layer: 10,
        wire: Wire::awg(28),
    };
    let tab = Table::build(&coil);
    let weight = MASS * G;
    let g0 = tab.zero_power_gap(weight);
    let (_, da) = tab.eval(&tab.a, g0);
    let (b, _) = tab.eval(&tab.b, g0);
    let (c, _) = tab.eval(&tab.c, g0);
    let (l, dl) = tab.eval(&tab.l, g0);
    let (_, dpsi) = tab.eval(&tab.psi, g0);
    let k = -da;
    let p = (k / MASS).sqrt();
    println!(
        "Corner: {MASS} kg on a 12 x 4 mm N42 disc; control coil {} turns AWG {} round it, R {:.3} ohm",
        coil.turns(),
        coil.wire.awg,
        tab.r
    );
    println!(
        "1. Linearized at the zero-power gap g0 = {:.3} mm (magnet alone carries {weight:.3} N)",
        g0 * 1e3
    );
    println!(
        "   negative stiffness k = -da/dg = {k:.1} N/m; open-loop pole p = sqrt(k/m) = {p:.1} 1/s \
         (doubling every {:.1} ms)",
        2f64.ln() / p * 1e3
    );
    println!(
        "   force per amp b = {b:.4} N/A, coil-coil c = {c:.2e} N/A^2; L = {:.3} mH, L/R = {:.3} ms",
        l * 1e3,
        l / tab.r * 1e3
    );
    println!(
        "   energy consistency: b vs -dPsi/dg {:.4} ({:+.1e}); c vs -L'/2 {:.2e} ({:+.1e})",
        -dpsi,
        b / -dpsi - 1.0,
        -0.5 * dl,
        c / (-0.5 * dl) - 1.0
    );

    println!("2. Open loop in sim-core: growth of a 10 um offset, measured vs predicted");
    let none = Gains {
        noise: 0.0,
        d_filter: 300.0,
        kp: 0.0,
        kd: 0.0,
        kpi: 0.0,
        kii: 0.0,
        zero_power: 0.0,
    };
    for dt in [5e-5, 2.5e-5] {
        let log = simulate(
            &tab,
            Drive::Open,
            none,
            g0 + 1e-5,
            g0,
            0.3,
            dt,
            (f64::MAX, 0.0),
        );
        let s = growth_rate(&log, g0);
        println!(
            "   coil open,    dt {:.0} us: {s:.2} 1/s vs p {p:.2} ({:+.1e}); ends at gap {:.2} mm",
            dt * 1e6,
            s / p - 1.0,
            log.gap.last().copied().unwrap_or(f64::NAN) * 1e3
        );
    }
    let ps = shorted_root(l, tab.r, b, k, MASS);
    for dt in [5e-5, 2.5e-5] {
        let log = simulate(
            &tab,
            Drive::Shorted,
            none,
            g0 + 1e-5,
            g0,
            0.3,
            dt,
            (f64::MAX, 0.0),
        );
        let s = growth_rate(&log, g0);
        println!(
            "   coil shorted, dt {:.0} us: {s:.2} 1/s vs cubic root {ps:.2} ({:+.1e})",
            dt * 1e6,
            s / ps - 1.0
        );
    }

    // Gains: gap loop at wn = 3 p (damping 0.7) through the force per amp; current loop at 1 kHz
    // by pole-zero cancellation of the coil's L/R.
    let wn = 3.0 * p;
    let wc = 2.0 * PI * 1000.0;
    let gains = Gains {
        noise: NOISE,
        d_filter: 100.0,
        kp: MASS * (wn * wn + p * p) / b,
        kd: 2.0 * 0.7 * wn * MASS / b,
        kpi: l * wc,
        kii: tab.r * wc,
        // Setpoint drift per amp: the zero-power loop settles in tau = b / (K k) = 0.5 s.
        zero_power: b / (k * 0.5),
    };
    println!(
        "3. Closed loop: gap PD at {SAMPLE:.0} Hz with {:.0} um sensor noise -> current command; PI \
         current loop -> coil voltage within +-{SUPPLY} V",
        NOISE * 1e6
    );
    println!(
        "   gains: Kp {:.1} A/m, Kd {:.3} A s/m (wn {wn:.0} 1/s), current loop {:.0} Hz",
        gains.kp,
        gains.kd,
        wc / (2.0 * PI)
    );
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out");
    std::fs::create_dir_all(&out).expect("create out/");
    for (label, zp, noise, file) in [
        ("hold, no noise", false, 0.0, "corner_hold_clean.svg"),
        ("gap hold", false, NOISE, "corner_hold.svg"),
        ("zero-power", true, NOISE, "corner_zero_power.svg"),
    ] {
        let gains = Gains { noise, ..gains };
        let dt = 2.5e-5;
        // Start 0.5 mm below g0, hold at g0; add 25 g payload at 1.5 s.
        let log = simulate(
            &tab,
            Drive::Control(zp),
            gains,
            g0 + 5e-4,
            g0,
            4.0,
            dt,
            (1.5, 0.025),
        );
        let window = |t0: f64, t1: f64, v: &[f64]| -> Vec<f64> {
            log.t
                .iter()
                .zip(v)
                .filter(|(t, _)| **t >= t0 && **t < t1)
                .map(|(_, x)| *x)
                .collect()
        };
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let rms = |v: &[f64]| {
            let m = mean(v);
            (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / v.len() as f64).sqrt()
        };
        let (g1, i1, p1) = (
            window(1.0, 1.5, &log.gap),
            window(1.0, 1.5, &log.current),
            window(1.0, 1.5, &log.power),
        );
        let (g2, i2, p2) = (
            window(3.5, 4.0, &log.gap),
            window(3.5, 4.0, &log.current),
            window(3.5, 4.0, &log.power),
        );
        println!(
            "   {label:<10}: before payload gap {:.3} mm (rms {:.1} um), current {:+.3} A, coil {:.1} mW; \
             after +25 g gap {:.3} mm, current {:+.3} A, coil {:.1} mW{}",
            mean(&g1) * 1e3,
            rms(&g1) * 1e6,
            mean(&i1),
            mean(&p1) * 1e3,
            mean(&g2) * 1e3,
            mean(&i2),
            mean(&p2) * 1e3,
            if log.t.last().copied().unwrap_or(0.0) < 3.99 {
                format!(
                    "  (lost at {:.4} s, gap {:.3} mm, current {:+.2} A)",
                    log.t.last().copied().unwrap_or(0.0),
                    log.gap.last().copied().unwrap_or(0.0) * 1e3,
                    log.current.last().copied().unwrap_or(0.0)
                )
            } else {
                String::new()
            }
        );
        if zp {
            let g_new = tab.zero_power_gap((MASS + 0.025) * G);
            println!(
                "     zero-power gap for the heavier corner, from the table: {:.3} mm",
                g_new * 1e3
            );
        }
        plot(&log, label, &out.join(file));
    }
    println!("   plots: {}", out.join("corner_*.svg").display());

    println!("4. Sensor noise vs derivative filter, zero-power control, 1.0-1.5 s (no payload)");
    println!("   noise  D filter  gap rms   current rms  coil power  held");
    for noise in [1e-6, 2e-6, 5e-6] {
        for d_filter in [100.0, 300.0] {
            let g = Gains {
                noise,
                d_filter,
                ..gains
            };
            let log = simulate(
                &tab,
                Drive::Control(true),
                g,
                g0 + 5e-4,
                g0,
                1.5,
                2.5e-5,
                (f64::MAX, 0.0),
            );
            let held = log.t.last().copied().unwrap_or(0.0) > 1.49;
            let sel: Vec<usize> = (0..log.t.len()).filter(|&k| log.t[k] >= 1.0).collect();
            let n = sel.len().max(1) as f64;
            let gm = sel.iter().map(|&k| log.gap[k]).sum::<f64>() / n;
            let grms = (sel.iter().map(|&k| (log.gap[k] - gm).powi(2)).sum::<f64>() / n).sqrt();
            let irms = (sel.iter().map(|&k| log.current[k].powi(2)).sum::<f64>() / n).sqrt();
            let pw = sel.iter().map(|&k| log.power[k]).sum::<f64>() / n;
            println!(
                "   {:>3.0} um  {d_filter:>4.0} Hz   {:>5.2} um  {irms:>7.3} A    {:>7.1} mW  {}",
                noise * 1e6,
                grms * 1e6,
                pw * 1e3,
                if held { "yes" } else { "no" }
            );
        }
    }
}

fn plot(log: &Log, title: &str, path: &std::path::Path) {
    let (w, h) = (900.0, 520.0);
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" font-family="sans-serif" font-size="13"><rect width="100%" height="100%" fill="#fff"/><text x="10" y="20" font-weight="bold">One corner, {title}: gap (mm) and coil current (A) vs time (s); +25 g payload at 1.5 s</text>"##
    );
    let t1 = log.t.last().copied().unwrap_or(1.0);
    let panel = |v: &[f64], y0: f64, ph: f64, color: &str, label: &str, svg: &mut String| {
        let (lo, hi) = v
            .iter()
            .fold((f64::MAX, f64::MIN), |a, &x| (a.0.min(x), a.1.max(x)));
        let pad = 0.05 * (hi - lo).max(1e-9);
        let (lo, hi) = (lo - pad, hi + pad);
        let x = |t: f64| 60.0 + (w - 80.0) * t / t1;
        let y = |val: f64| y0 + ph - ph * (val - lo) / (hi - lo);
        let pts: Vec<String> = log
            .t
            .iter()
            .zip(v)
            .step_by(4)
            .map(|(&t, &val)| format!("{:.1},{:.1}", x(t), y(val)))
            .collect();
        let _ = writeln!(
            svg,
            r##"<rect x="60" y="{y0}" width="{}" height="{ph}" fill="none" stroke="#ccc"/>"##,
            w - 80.0
        );
        let _ = writeln!(
            svg,
            r#"<polyline points="{}" fill="none" stroke="{color}" stroke-width="1.2"/>"#,
            pts.join(" ")
        );
        let _ = writeln!(
            svg,
            r##"<text x="62" y="{:.0}" fill="#444">{label}: {hi:.4} .. {lo:.4}</text>"##,
            y0 + 14.0
        );
    };
    let gap_mm: Vec<f64> = log.gap.iter().map(|g| g * 1e3).collect();
    panel(&gap_mm, 35.0, 220.0, "#2e86c1", "gap mm", &mut svg);
    panel(&log.current, 280.0, 220.0, "#c0392b", "current A", &mut svg);
    let _ = writeln!(svg, "</svg>");
    std::fs::write(path, svg).expect("write SVG");
}
