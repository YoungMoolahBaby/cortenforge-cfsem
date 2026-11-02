//! `cart`: the whole demonstrator cart in sim-core, four hybrid corners under the rail.
//!
//! One body on six joints (slide along the track, sideways, heave, yaw, pitch, roll; sideways and
//! yaw held by the guidance coils of `guide`). Sites mark the four magnet faces. Inside
//! `cb_passive` each corner's gap and gap rate come from its site and Jacobian, and the cfsem
//! force table's pull F(g, I) is added to the joint forces as J^T F; the linear motor's and the
//! guidance coils' wrenches, a sideways shove and an off-centre payload go in the same way. Coil
//! currents and the controllers run each step outside, sharing state with the callback.
//!
//! Control is modal: heave, pitch and roll each get a gap PD plus a zero-power loop; the fourth
//! corner pattern (warp) cannot move a rigid body, so its current is held at zero.

use std::f64::consts::PI;
use std::fmt::Write as _;
use std::sync::{Arc, Mutex};

use cortenforge::sim::core::types::callbacks::Callback;
use cortenforge::sim::core::{Model, mj_jac, mj_jac_site};
use cortenforge::sim::mjcf::load_model;
use nalgebra::Vector3;

use crate::corner::{Table, corner_coil};
use crate::guide::GuideTable;
use crate::motor::MotorTable;

const G: f64 = 9.81;
const MASS: f64 = 0.5; // kg, cart with its nominal payload
const HALF: (f64, f64, f64) = (0.07, 0.04, 0.015); // m, body box half-extents
const LX: f64 = 0.05;
const LY: f64 = 0.025;
/// Corners (x, y): 0 rear-right, 1 front-right, 2 front-left, 3 rear-left.
const CORNERS: [(f64, f64); 4] = [(-LX, -LY), (LX, -LY), (LX, LY), (-LX, LY)];
/// Mode shapes over the corners: heave, pitch (front up), roll (left up). The fourth pattern,
/// warp (+1, -1, +1, -1), cannot move a rigid body; corner currents are built from these three
/// only, so the warp current is zero by construction.
const SH: [f64; 4] = [1.0, 1.0, 1.0, 1.0];
const SP: [f64; 4] = [-1.0, 1.0, 1.0, -1.0];
const SR: [f64; 4] = [-1.0, -1.0, 1.0, 1.0];
const SUPPLY: f64 = 12.0;
const I_MAX: f64 = 10.0;
const NOISE: f64 = 1e-6;
const SAMPLE: f64 = 2000.0;
const D_FILTER: f64 = 100.0;
const DT: f64 = 2.5e-5;
const RAMP: f64 = 0.3;
const START_GAP: f64 = 0.0065;
/// Extra load set down at t_payload, at (x, y) on the cart top.
const PAYLOAD: (f64, f64, f64) = (0.1, 0.03, 0.015);
const T_PAYLOAD: f64 = 1.0;
/// The trip: from x = 0 to the station, trapezoidal speed profile.
const STATION: f64 = 1.2;
const V_CRUISE: f64 = 0.3;
const ACCEL: f64 = 0.5;
const T_DEPART: f64 = 2.5;
/// Sideways start: offset (m) and yaw (rad); a 0.3 N shove at the front edge for 50 ms at T_PUSH.
const Y0: f64 = 0.001;
const YAW0: f64 = 0.5 * PI / 180.0;
const T_PUSH: f64 = 5.0;
const PUSH: f64 = 0.3;
/// Sideways sensors at x = +-SENSOR_X (m), noise 5 um; lip clearance 2 mm at the body corners.
const SENSOR_X: f64 = 0.06;
const LAT_NOISE: f64 = 1e-6;
/// Guidance current limit (A) and loop bandwidth (Hz). A first try at 3 A and 4 Hz let the 0.3 N
/// shove reach the lips: 3 A in both blocks makes at most about 18 mN m of yaw, under the shove's 21.
const I_GUIDE: f64 = 8.0;
const F_GUIDE: f64 = 8.0;
const LIP: f64 = 0.002;

/// State shared between the physics callback and the controller.
#[derive(Default)]
struct Shared {
    current: [f64; 4],
    gap: [f64; 4],
    gap_rate: [f64; 4],
    /// (A) linear-motor thrust current.
    iq: f64,
    payload_on: bool,
    /// (A) guidance currents of the front and rear blocks.
    i_front: f64,
    i_rear: f64,
    /// (N) a sideways shove at the cart's front edge (a disturbance).
    push: f64,
}

/// Reference position and speed along the track at time t.
fn trip(t: f64) -> (f64, f64, f64) {
    let t = t - T_DEPART;
    if t <= 0.0 {
        return (0.0, 0.0, 0.0);
    }
    let ta = V_CRUISE / ACCEL;
    let d_acc = 0.5 * ACCEL * ta * ta;
    let t_cruise = (STATION - 2.0 * d_acc) / V_CRUISE;
    if t < ta {
        (0.5 * ACCEL * t * t, ACCEL * t, ACCEL)
    } else if t < ta + t_cruise {
        (d_acc + V_CRUISE * (t - ta), V_CRUISE, 0.0)
    } else if t < 2.0 * ta + t_cruise {
        let u = t - ta - t_cruise;
        (
            d_acc + V_CRUISE * t_cruise + V_CRUISE * u - 0.5 * ACCEL * u * u,
            V_CRUISE - ACCEL * u,
            -ACCEL,
        )
    } else {
        (STATION, 0.0, 0.0)
    }
}

/// Motor wrench per amp at cart position x, interpolated linearly between tables at their gaps.
fn motor_wrench(motor: &[(f64, MotorTable)], x: f64, gap: f64) -> [f64; 4] {
    let k = motor
        .iter()
        .position(|(g, _)| *g >= gap)
        .unwrap_or(motor.len() - 1)
        .clamp(1, motor.len() - 1);
    let ((g0, t0), (g1, t1)) = (&motor[k - 1], &motor[k]);
    let u = ((gap - g0) / (g1 - g0)).clamp(0.0, 1.0);
    let (a, b) = (t0.at(x), t1.at(x));
    [0, 1, 2, 3].map(|c| a[c] * (1.0 - u) + b[c] * u)
}

fn build_model(
    tab: Arc<Table>,
    shared: Arc<Mutex<Shared>>,
    motor: Arc<Vec<(f64, MotorTable)>>,
    guide: Arc<GuideTable>,
) -> Model {
    let sites: String = CORNERS
        .iter()
        .enumerate()
        .map(|(k, (x, y))| format!(r#"<site name="c{k}" pos="{x} {y} {}"/>"#, HALF.2))
        .collect();
    let xml = format!(
        r#"<mujoco model="cart">
  <option timestep="{DT}" gravity="0 0 -{G}"/>
  <worldbody>
    <body name="cart" pos="0 0 {}">
      <joint name="x" type="slide" axis="1 0 0"/>
      <joint name="y" type="slide" axis="0 1 0"/>
      <joint name="z" type="slide" axis="0 0 1"/>
      <joint name="yaw" type="hinge" axis="0 0 1"/>
      <joint name="pitch" type="hinge" axis="0 -1 0"/>
      <joint name="roll" type="hinge" axis="1 0 0"/>
      <geom type="box" size="{} {} {}" mass="{MASS}"/>
      {sites}
    </body>
  </worldbody>
</mujoco>"#,
        -START_GAP - HALF.2,
        HALF.0,
        HALF.1,
        HALF.2
    );
    let mut model = load_model(&xml).expect("cart model loads");
    model.cb_passive = Some(Callback(Arc::new(
        move |m: &Model, d: &mut cortenforge::sim::core::Data| {
            let mut s = shared.lock().expect("shared state");
            let body = 1;
            for k in 0..4 {
                let (jp, _) = mj_jac_site(m, d, k);
                let gap = -d.site_xpos[k].z;
                let vz: f64 = (0..m.nv).map(|j| jp[(2, j)] * d.qvel[j]).sum();
                s.gap[k] = gap;
                s.gap_rate[k] = -vz;
                let f = tab.force(gap, s.current[k]);
                for j in 0..m.nv {
                    d.qfrc_passive[j] += jp[(2, j)] * f;
                }
            }
            // The linear motor: its wrench per amp (thrust, lift, moments about x and y, at the
            // centre of mass) from the cfsem tables, at this position and mean gap.
            let mean_gap = s.gap.iter().sum::<f64>() / 4.0;
            let w = motor_wrench(&motor, d.qpos[0], mean_gap);
            let com = d.xipos[body];
            let (jc, jr) = mj_jac(m, d, body, &com);
            // The payload's weight at its point.
            let at = |p: Vector3<f64>| d.xpos[body] + d.xmat[body] * p;
            let payload = if s.payload_on { PAYLOAD.0 * G } else { 0.0 };
            let pp = at(Vector3::new(PAYLOAD.1, PAYLOAD.2, HALF.2));
            let (jw, _) = mj_jac(m, d, body, &pp);
            // Guidance: the front and rear blocks' wrenches per amp, from the cfsem tables.
            let (gf, gr) = guide.at(d.qpos[0]);
            let g: [f64; 6] = [0, 1, 2, 3, 4, 5].map(|c| s.i_front * gf[c] + s.i_rear * gr[c]);
            // A sideways shove at the front edge.
            let pf = at(Vector3::new(HALF.0, 0.0, 0.0));
            let (jf, _) = mj_jac(m, d, body, &pf);
            for j in 0..m.nv {
                d.qfrc_passive[j] += s.iq
                    * (jc[(0, j)] * w[0]
                        + jc[(2, j)] * w[1]
                        + jr[(0, j)] * w[2]
                        + jr[(1, j)] * w[3])
                    - jw[(2, j)] * payload
                    + jc[(0, j)] * g[0]
                    + jc[(1, j)] * g[1]
                    + jc[(2, j)] * g[2]
                    + jr[(0, j)] * g[3]
                    + jr[(1, j)] * g[4]
                    + jr[(2, j)] * g[5]
                    + jf[(1, j)] * s.push;
            }
        },
    )));
    model
}

/// Per-mode control: PD gains and zero-power drift, for an equivalent per-corner mass.
#[derive(Clone, Copy)]
struct ModeGains {
    kp: f64,
    kd: f64,
    zp: f64,
}

fn mode_gains(m_eq: f64, k: f64, b: f64) -> ModeGains {
    let p = (k / m_eq).sqrt();
    let wn = 3.0 * p;
    ModeGains {
        kp: m_eq * (wn * wn + p * p) / b,
        kd: 2.0 * 0.7 * wn * m_eq / b,
        zp: b / (k * 0.5),
    }
}

struct Log {
    t: Vec<f64>,
    y: Vec<f64>,
    yaw: Vec<f64>,
    guide: Vec<[f64; 2]>,
    iq: Vec<f64>,
    stator: Vec<f64>,
    x: Vec<f64>,
    gap: Vec<[f64; 4]>,
    current: Vec<[f64; 4]>,
    power: Vec<f64>,
}

/// The zero-power plane: heave, pitch, roll mode gaps at which the four magnets alone hold the
/// cart and payload in force and moment balance (Newton on the force table).
fn zero_power_plane(tab: &Table, payload: bool) -> [f64; 3] {
    let w = MASS * G + if payload { PAYLOAD.0 * G } else { 0.0 };
    let (mx, my) = if payload {
        (PAYLOAD.1 * PAYLOAD.0 * G, PAYLOAD.2 * PAYLOAD.0 * G)
    } else {
        (0.0, 0.0)
    };
    let resid = |v: [f64; 3]| -> [f64; 3] {
        let mut r = [-w, -mx, -my];
        for k in 0..4 {
            let g = v[0] + SP[k] * v[1] + SR[k] * v[2];
            let f = tab.eval(&tab.a, g).0;
            r[0] += f;
            r[1] += CORNERS[k].0 * f;
            r[2] += CORNERS[k].1 * f;
        }
        // The centre of mass sits HALF.2 below the magnet faces. Tilting front-up by
        // alpha = -P / LX (roll left-up by beta = -R / LY) carries the faces back (right) by
        // HALF.2 alpha relative to it, so the cart's weight needs HALF.2 alpha W more front
        // (left) moment from the magnets.
        let wc = MASS * G;
        r[1] -= HALF.2 * (-v[1] / LX) * wc;
        r[2] -= HALF.2 * (-v[2] / LY) * wc;
        r
    };
    let mut v = [tab.zero_power_gap(w / 4.0), 0.0, 0.0];
    for _ in 0..30 {
        let r0 = resid(v);
        let h = 1e-7;
        let mut jm = [[0.0; 3]; 3];
        for c in 0..3 {
            let mut vp = v;
            vp[c] += h;
            let rp = resid(vp);
            for r in 0..3 {
                jm[r][c] = (rp[r] - r0[r]) / h;
            }
        }
        let d = nalgebra::Matrix3::from_fn(|r, c| jm[r][c])
            .lu()
            .solve(&nalgebra::Vector3::new(-r0[0], -r0[1], -r0[2]))
            .expect("solvable");
        for c in 0..3 {
            v[c] += d[c];
        }
    }
    v
}

pub fn run() {
    let coil = corner_coil();
    let tab = Arc::new(Table::build(&coil));
    let shared = Arc::new(Mutex::new(Shared::default()));
    // Motor tables at four gaps (the cart rides between about 5 and 6.5 mm).
    let motor: Arc<Vec<(f64, MotorTable)>> = Arc::new(
        [0.005, 0.0055, 0.006, 0.0065]
            .iter()
            .map(|&g| (g, MotorTable::build(g, 24)))
            .collect(),
    );
    let kf = motor_wrench(&motor, 0.0, 0.005964)[0];
    let copper = motor[2].1.copper;
    // Guidance tables at the nominal gap (the force changes about 1.7% a mm off-centre).
    let guide = Arc::new(GuideTable::build(0.005964, 24));
    let model = Arc::new(build_model(
        tab.clone(),
        shared.clone(),
        motor.clone(),
        guide.clone(),
    ));
    let mut data = model.make_data();
    let inertia = model.body_inertia[1];
    // Start 1 mm off the track centre and turned 0.5 deg.
    data.qpos[1] = Y0;
    data.qpos[3] = YAW0;
    // Sideways and yaw control: PD at F_GUIDE Hz, damping 0.8, on the cart's mass and yaw inertia.
    let wl = 2.0 * PI * F_GUIDE;
    let (kpy, kdy) = (MASS * wl * wl, 2.0 * 0.8 * wl * MASS);
    let (kpz, kdz) = (inertia.z * wl * wl, 2.0 * 0.8 * wl * inertia.z);
    let (mut lat_last, mut lat_rate) = ([Y0, YAW0], [0.0; 2]);
    let mut lat_cmd = [0.0; 2];
    let g0 = tab.zero_power_gap(MASS * G / 4.0);
    let k = -tab.eval(&tab.a, g0).1;
    let b = tab.eval(&tab.b, g0).0;
    // Equivalent per-corner masses: heave M/4, pitch I_yy / (4 lx^2), roll I_xx / (4 ly^2).
    let m_h = MASS / 4.0;
    let m_p = inertia.y / (4.0 * LX * LX);
    let m_r = inertia.x / (4.0 * LY * LY);
    let gh = mode_gains(m_h, k, b);
    let gp = mode_gains(m_p, k, b);
    let gr = mode_gains(m_r, k, b);
    println!(
        "Cart: {MASS} kg, body {:.0} x {:.0} x {:.0} mm (I_xx {:.2e}, I_yy {:.2e} kg m^2), corners at +-{:.0} / +-{:.0} mm",
        2e3 * HALF.0,
        2e3 * HALF.1,
        2e3 * HALF.2,
        inertia.x,
        inertia.y,
        LX * 1e3,
        LY * 1e3
    );
    println!(
        "   per-corner equivalent mass: heave {m_h:.4}, pitch {m_p:.4}, roll {m_r:.4} kg; open-loop poles {:.1} / {:.1} / {:.1} 1/s",
        (k / m_h).sqrt(),
        (k / m_p).sqrt(),
        (k / m_r).sqrt()
    );

    // Controller state.
    let mut i_int = [0.0; 4];
    let mut i_cmd = [0.0; 4];
    let mut refs = [g0, 0.0, 0.0]; // heave, pitch, roll mode setpoints
    let mut last = [START_GAP, 0.0, 0.0];
    let mut rate = [0.0; 3];
    let mut mode_cmd = [0.0; 3];
    let steps_per_sample = ((1.0 / SAMPLE) / DT).round() as usize;
    let mut rng = 0x9e37_79b9_7f4a_7c15_u64;
    let mut gauss = || {
        let mut u = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            (rng >> 11) as f64 / (1u64 << 53) as f64
        };
        let (a, c) = (u().max(1e-300), u());
        (-2.0 * a.ln()).sqrt() * (2.0 * PI * c).cos()
    };
    let mut log = Log {
        t: vec![],
        y: vec![],
        yaw: vec![],
        guide: vec![],
        iq: vec![],
        stator: vec![],
        x: vec![],
        gap: vec![],
        current: vec![],
        power: vec![],
    };
    let t_end = 9.0;
    let n = (t_end / DT).round() as usize;
    let mut lost = None;
    for step in 0..n {
        let t = step as f64 * DT;
        // Physics first (the callback refreshes gaps), then the controller acts on them.
        data.step(&model).expect("step");
        let (gap, gap_rate) = {
            let s = shared.lock().expect("shared state");
            (s.gap, s.gap_rate)
        };
        if gap.iter().any(|g| !(0.0015..0.02).contains(g)) {
            lost = Some((t, gap));
            break;
        }
        if step % steps_per_sample == 0 {
            let ts = 1.0 / SAMPLE;
            let meas: Vec<f64> = gap.iter().map(|g| g + NOISE * gauss()).collect();
            let modes = [
                (0..4).map(|c| SH[c] * meas[c]).sum::<f64>() / 4.0,
                (0..4).map(|c| SP[c] * meas[c]).sum::<f64>() / 4.0,
                (0..4).map(|c| SR[c] * meas[c]).sum::<f64>() / 4.0,
            ];
            let alpha = 1.0 - (-2.0 * PI * D_FILTER * ts).exp();
            for (mi, gains) in [gh, gp, gr].iter().enumerate() {
                rate[mi] += alpha * ((modes[mi] - last[mi]) / ts - rate[mi]);
                last[mi] = modes[mi];
                if t > RAMP {
                    refs[mi] -= gains.zp * mode_cmd[mi] * ts;
                }
                let target = if mi == 0 && t < RAMP {
                    START_GAP + (refs[0] - START_GAP) * t / RAMP
                } else {
                    refs[mi]
                };
                mode_cmd[mi] = gains.kp * (modes[mi] - target) + gains.kd * rate[mi];
            }
            for c in 0..4 {
                i_cmd[c] =
                    (mode_cmd[0] + SP[c] * mode_cmd[1] + SR[c] * mode_cmd[2]).clamp(-I_MAX, I_MAX);
            }
            // Sideways: two sensors at x = +-SENSOR_X read the cart's sideways position there.
            let (y, yaw) = (data.qpos[1], data.qpos[3]);
            let yf = y + SENSOR_X * yaw.sin() + LAT_NOISE * gauss();
            let yr = y - SENSOR_X * yaw.sin() + LAT_NOISE * gauss();
            let lat = [0.5 * (yf + yr), ((yf - yr) / (2.0 * SENSOR_X)).asin()];
            let alpha = 1.0 - (-2.0 * PI * 50.0 * ts).exp();
            for c in 0..2 {
                lat_rate[c] += alpha * ((lat[c] - lat_last[c]) / ts - lat_rate[c]);
            }
            lat_last = lat;
            lat_cmd = [
                -(kpy * lat[0] + kdy * lat_rate[0]),
                -(kpz * lat[1] + kdz * lat_rate[1]),
            ];
        }
        // Current loops and coil equations, per corner.
        let mut s = shared.lock().expect("shared state");
        let mut power = 0.0;
        for c in 0..4 {
            let i = s.current[c];
            let e = i_cmd[c] - i;
            let v_raw = 1.22e-4 * 2.0 * PI * 1000.0 * e + tab.r * 2.0 * PI * 1000.0 * i_int[c];
            if v_raw.abs() < SUPPLY || v_raw.signum() != e.signum() {
                i_int[c] += e * DT;
            }
            let v = v_raw.clamp(-SUPPLY, SUPPLY);
            let (l, dl) = tab.eval(&tab.l, gap[c]);
            let (_, dpsi) = tab.eval(&tab.psi, gap[c]);
            s.current[c] = i + (v - tab.r * i - (i * dl + dpsi) * gap_rate[c]) / l * DT;
            power += s.current[c] * s.current[c] * tab.r;
        }
        // Propulsion: track the trip profile (mass feedforward plus PD), and the payload.
        let (xr, vr, ar) = trip(t);
        let f_cmd = MASS * ar + 40.0 * (xr - data.qpos[0]) + 12.0 * (vr - data.qvel[0]);
        s.iq = (f_cmd / kf).clamp(-3.0, 3.0);
        s.payload_on = t >= T_PAYLOAD;
        // Guidance currents: solve the front and rear blocks' sideways force and yaw moment per
        // amp (from the tables at this position) for the commanded force and moment.
        let (gf, gr) = guide.at(data.qpos[0]);
        let det = gf[1] * gr[5] - gr[1] * gf[5];
        s.i_front = ((lat_cmd[0] * gr[5] - gr[1] * lat_cmd[1]) / det).clamp(-I_GUIDE, I_GUIDE);
        s.i_rear = ((gf[1] * lat_cmd[1] - lat_cmd[0] * gf[5]) / det).clamp(-I_GUIDE, I_GUIDE);
        s.push = if (T_PUSH..T_PUSH + 0.05).contains(&t) {
            PUSH
        } else {
            0.0
        };
        if step % 40 == 0 {
            log.t.push(t);
            log.y.push(data.qpos[1]);
            log.yaw.push(data.qpos[3]);
            log.guide.push([s.i_front, s.i_rear]);
            log.iq.push(s.iq);
            log.stator.push(s.iq * s.iq * copper);
            log.x.push(data.qpos[0]);
            log.gap.push(gap);
            log.current.push(s.current);
            log.power.push(power);
        }
    }
    if let Some((t, gap)) = lost {
        println!(
            "   LOST at {t:.3} s, gaps {:?} mm",
            gap.map(|g| (g * 1e4).round() / 10.0)
        );
    }

    // Phases.
    let phase = |t0: f64, t1: f64| -> Vec<usize> {
        (0..log.t.len())
            .filter(|&i| log.t[i] >= t0 && log.t[i] < t1)
            .collect()
    };
    let mean_gap = |idx: &[usize], c: usize| {
        idx.iter().map(|&i| log.gap[i][c]).sum::<f64>() / idx.len() as f64
    };
    let mean_i = |idx: &[usize], c: usize| {
        idx.iter().map(|&i| log.current[i][c]).sum::<f64>() / idx.len() as f64
    };
    let mean_p = |idx: &[usize]| idx.iter().map(|&i| log.power[i]).sum::<f64>() / idx.len() as f64;
    println!(
        "Timeline: lift-off ramp 0-{RAMP} s; +{:.0} g set down at ({:+.0}, {:+.0}) mm at {T_PAYLOAD} s; depart at {T_DEPART} s for the station {STATION} m away ({V_CRUISE} m/s, {ACCEL} m/s^2)",
        PAYLOAD.0 * 1e3,
        PAYLOAD.1 * 1e3,
        PAYLOAD.2 * 1e3
    );
    let show = |name: &str, idx: &[usize]| {
        println!(
            "   {name:<28} gaps {:.3} / {:.3} / {:.3} / {:.3} mm, currents {:+.3} / {:+.3} / {:+.3} / {:+.3} A, coils {:.0} mW",
            mean_gap(idx, 0) * 1e3,
            mean_gap(idx, 1) * 1e3,
            mean_gap(idx, 2) * 1e3,
            mean_gap(idx, 3) * 1e3,
            mean_i(idx, 0),
            mean_i(idx, 1),
            mean_i(idx, 2),
            mean_i(idx, 3),
            mean_p(idx) * 1e3
        );
    };
    println!("   (corners: rear-right / front-right / front-left / rear-left)");
    show("hovering, 0.7-1.0 s", &phase(0.7, 1.0));
    show("with payload, 2.0-2.5 s", &phase(2.0, 2.5));
    let trip_idx = phase(T_DEPART, 8.0);
    let worst = trip_idx
        .iter()
        .map(|&i| {
            (0..4)
                .map(|c| (log.gap[i][c] - mean_gap(&phase(2.0, 2.5), c)).abs())
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);
    println!(
        "   during the trip ({T_DEPART}-8.0 s): largest corner gap excursion {:.1} um, mean coil power {:.0} mW",
        worst * 1e6,
        mean_p(&trip_idx) * 1e3
    );
    show("at the station, 8.5-9.0 s", &phase(8.5, 9.0));
    let xf = log.x.last().copied().unwrap_or(f64::NAN);
    println!(
        "   stopped at x = {:.4} m (station {STATION} m, {:+.1} mm)",
        xf,
        (xf - STATION) * 1e3
    );
    let dt_log = log.t.get(1).copied().unwrap_or(0.0) - log.t.first().copied().unwrap_or(0.0);
    let peak = log.iq.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    let energy: f64 = log.stator.iter().sum::<f64>() * dt_log;
    println!(
        "   linear motor: thrust constant {kf:.4} N/A, copper {copper:.2} W/A^2; peak thrust current {peak:.3} A, \
         stator energy for the trip {energy:.3} J (peak {:.2} W)",
        log.stator.iter().fold(0.0_f64, |a, &v| a.max(v))
    );
    // Sideways: the body corners (x = +-HALF.0) against the lips LIP away.
    let corner = |i: usize| log.y[i].abs() + HALF.0 * log.yaw[i].sin().abs();
    let worst = |idx: &[usize]| {
        idx.iter()
            .map(|&i| (corner(i), log.t[i]))
            .fold((0.0, 0.0), |a, b| if b.0 > a.0 { b } else { a })
    };
    let all: Vec<usize> = (0..log.t.len()).collect();
    let (w_all, t_all) = worst(&all);
    let (w_push, _) = worst(&phase(T_PUSH, T_PUSH + 1.0));
    let settled = phase(2.0, 2.5);
    let ms = |v: &dyn Fn(usize) -> f64, idx: &[usize]| {
        (idx.iter().map(|&i| v(i) * v(i)).sum::<f64>() / idx.len() as f64).sqrt()
    };
    let gcur = |i: usize| log.guide[i][0].hypot(log.guide[i][1]);
    let g_energy: f64 = log
        .guide
        .iter()
        .map(|g| (g[0] * g[0] + g[1] * g[1]) * guide.copper)
        .sum::<f64>()
        * dt_log;
    println!(
        "   sideways: start {:.1} mm off and {:.1} deg; settled by 2.0-2.5 s to {:.1} um rms, {:.4} deg rms; \
         worst body-corner excursion {:.3} mm at {:.2} s; after the {PUSH} N x 50 ms shove {:.3} mm; \
         lips at {:.1} mm: {}",
        Y0 * 1e3,
        YAW0.to_degrees(),
        ms(&|i| log.y[i], &settled) * 1e6,
        ms(&|i| log.yaw[i], &settled).to_degrees(),
        w_all * 1e3,
        t_all,
        w_push * 1e3,
        LIP * 1e3,
        if w_all < LIP { "no contact" } else { "CONTACT" }
    );
    println!(
        "   guidance: peak {:.2} A, rms {:.3} A over the run, {:.3} A while hovering (0.7-1.0 s, {:.0} mW: sensor noise); \n         copper {:.2} W/A^2 a block, {:.3} J for the whole run",
        log.guide
            .iter()
            .fold(0.0_f64, |a, g| a.max(g[0].abs()).max(g[1].abs())),
        ms(&gcur, &all),
        ms(&gcur, &phase(0.7, 1.0)),
        ms(&gcur, &phase(0.7, 1.0)).powi(2) * guide.copper * 1e3,
        guide.copper,
        g_energy
    );
    for (label, payload, idx) in [
        ("no payload", false, phase(0.7, 1.0)),
        ("with payload", true, phase(8.5, 9.0)),
    ] {
        let v = zero_power_plane(&tab, payload);
        let pred: Vec<f64> = (0..4).map(|c| v[0] + SP[c] * v[1] + SR[c] * v[2]).collect();
        let worst = (0..4)
            .map(|c| (mean_gap(&idx, c) - pred[c]).abs())
            .fold(0.0, f64::max);
        println!(
            "   zero-power plane from the force table, {label}: {:.3} / {:.3} / {:.3} / {:.3} mm; sim within {:.1} um",
            pred[0] * 1e3,
            pred[1] * 1e3,
            pred[2] * 1e3,
            pred[3] * 1e3,
            worst * 1e6
        );
    }
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out");
    std::fs::create_dir_all(&out).expect("create out/");
    plot(&log, &out.join("cart.svg"));
    println!("   plot: {}", out.join("cart.svg").display());
}

fn plot(log: &Log, path: &std::path::Path) {
    let (w, h) = (1000.0, 975.0);
    let t1 = log.t.last().copied().unwrap_or(1.0);
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" font-family="sans-serif" font-size="13"><rect width="100%" height="100%" fill="#fff"/><text x="10" y="20" font-weight="bold">Cart: position (m), corner gaps (mm), corner currents (A), sideways (mm) vs time (s)</text>"##
    );
    let colors = ["#c0392b", "#e67e22", "#27ae60", "#2e86c1"];
    let mut panel = |series: Vec<Vec<f64>>, y0: f64, ph: f64, label: &str| {
        let (lo, hi) = series
            .iter()
            .flatten()
            .fold((f64::MAX, f64::MIN), |a, &x| (a.0.min(x), a.1.max(x)));
        let pad = 0.05 * (hi - lo).max(1e-9);
        let (lo, hi) = (lo - pad, hi + pad);
        let x = |t: f64| 60.0 + (w - 80.0) * t / t1;
        let y = |v: f64| y0 + ph - ph * (v - lo) / (hi - lo);
        let _ = writeln!(
            svg,
            r##"<rect x="60" y="{y0}" width="{}" height="{ph}" fill="none" stroke="#ccc"/><text x="62" y="{:.0}" fill="#444">{label}: {lo:.4} .. {hi:.4}</text>"##,
            w - 80.0,
            y0 + 14.0
        );
        for (k, s) in series.iter().enumerate() {
            let pts: Vec<String> = log
                .t
                .iter()
                .zip(s)
                .map(|(&t, &v)| format!("{:.1},{:.1}", x(t), y(v)))
                .collect();
            let _ = writeln!(
                svg,
                r#"<polyline points="{}" fill="none" stroke="{}" stroke-width="1"/>"#,
                pts.join(" "),
                colors[k % 4]
            );
        }
    };
    panel(vec![log.x.clone()], 35.0, 160.0, "position m");
    panel(
        (0..4)
            .map(|c| log.gap.iter().map(|g| g[c] * 1e3).collect())
            .collect(),
        215.0,
        260.0,
        "corner gaps mm",
    );
    panel(
        (0..4)
            .map(|c| log.current.iter().map(|i| i[c]).collect())
            .collect(),
        495.0,
        240.0,
        "corner currents A",
    );
    // Sideways: the centre's offset and the worst body corner's, against the lips.
    panel(
        vec![
            log.y.iter().map(|y| y * 1e3).collect(),
            log.y
                .iter()
                .zip(&log.yaw)
                .map(|(y, a)| (y.abs() + HALF.0 * a.sin().abs()) * 1e3)
                .collect(),
            vec![LIP * 1e3; log.t.len()],
        ],
        755.0,
        200.0,
        "sideways mm (centre, worst corner, lip)",
    );
    let _ = writeln!(svg, "</svg>");
    std::fs::write(path, svg).expect("write SVG");
}
