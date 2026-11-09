//! `swim`: the desk swimmer in glycerol. How far it goes a turn and how hard it is to turn, the
//! field rate where the magnet can no longer keep up (step-out) and the slipping beyond it,
//! whether gravity wins, and the Reynolds numbers of the desk demo and a cell-sized swimmer.
//!
//! Below step-out the swimmer turns with the field, its magnet lagging by phi with
//! m B sin(phi) = D omega (D the turning drag). Above it the lag runs round (Adler's equation,
//! d phi / dt = omega - omega_so sin phi): the swimmer turns on average at
//! omega - sqrt(omega^2 - omega_so^2), and Stokes flow being linear, its speed follows.

use std::f64::consts::PI;

use nalgebra::{Matrix3, Vector3};

use crate::check::{GLYCEROL_MU, GLYCEROL_RHO};
use crate::rft::{Swim, advance, swim, swim_locked};
use crate::swimmer::Swimmer;

const G: f64 = 9.81;
/// Water at 20 C, for the cell-sized comparison: viscosity (Pa s) and density (kg/m^3).
const WATER_MU: f64 = 1.0e-3;
const WATER_RHO: f64 = 998.0;

/// Average spin (Hz) of a swimmer whose step-out is `f_so`, in a field turning at `f`.
pub fn spin(f: f64, f_so: f64) -> f64 {
    if f <= f_so {
        f
    } else {
        f - (f * f - f_so * f_so).sqrt()
    }
}

pub fn run() {
    let s = Swimmer::desk();
    let h = s.helix;
    let mass = s.mass();
    let m = s.moment();
    println!("swim: the desk swimmer in glycerol ({GLYCEROL_MU} Pa s, {GLYCEROL_RHO} kg/m^3)");
    println!(
        "   helix: R {:.2} mm, pitch {:.1} mm ({:.1} deg), wire {:.1} mm across, {} turns",
        h.radius * 1e3,
        h.pitch * 1e3,
        h.pitch_angle_deg(),
        2.0 * h.wire * 1e3,
        h.turns
    );
    println!(
        "   head: {:.1} mm across, {:.1} mm long; magnet {:.0} x {:.0} mm across its axis, Br {} T, m = {:.3e} A m^2",
        2.0 * s.head_radius * 1e3,
        (s.head_len + 2.0 * s.head_radius) * 1e3,
        2.0 * s.mag_radius * 1e3,
        s.mag_len * 1e3,
        s.br,
        m
    );
    println!(
        "   length {:.2} mm, volume {:.1} mm^3, mass {:.1} mg (resin {}, magnet {} kg/m^3)",
        s.length() * 1e3,
        mass.volume * 1e9,
        mass.mass * 1e6,
        s.rho_body,
        s.rho_magnet
    );
    let model = s.model(GLYCEROL_MU);
    println!(
        "   model: {} points at 16 round, extrapolated from 12 and 16 ({:.0} s)",
        model.points, model.seconds
    );

    // -----------------------------------------------------------------------------------------
    println!("1. Turning it, per turn of the swimmer");
    let per_hz = |w: &Swim| (w.u.z * 2.0 * PI, w.torque * 2.0 * PI);
    let row = |what: &str, w: &Swim| {
        let (v, t) = per_hz(w);
        println!(
            "   {what:<34} advance {:.4} of the pitch, {:.3} mm/s per Hz, torque {:.3e} N m per Hz, wobble {:.3}",
            advance(&h, w),
            v * 1e3,
            t,
            (w.wobble.0.powi(2) + w.wobble.1.powi(2)).sqrt()
        );
    };
    let held = swim_locked(&model.r);
    let free = swim(&model.r);
    row("axis held by the field", &held);
    row("axis free to wobble", &free);
    row("axis held, 12 round alone", &swim_locked(&model.r12));
    row("axis held, 16 round alone", &swim_locked(&model.r16));
    println!(
        "   (the real swimmer lies between held and free: the field also resists tilting; held is used below)"
    );

    // -----------------------------------------------------------------------------------------
    println!("2. Step-out: the fastest the magnet can turn it, f = m B / (2 pi D)");
    let (v_hz, t_hz) = per_hz(&held);
    let d = held.torque; // N m per rad/s
    println!(
        "   {:>6} {:>12} {:>14} {:>14} {:>16}",
        "B (mT)", "step-out Hz", "top mm/s", "lengths / s", "free axis: Hz"
    );
    for b in [1.0, 2.0, 5.0, 10.0] {
        let f_so = m * b * 1e-3 / (2.0 * PI * d);
        let f_free = m * b * 1e-3 / (2.0 * PI * free.torque);
        println!(
            "   {b:>6.0} {f_so:>12.2} {:>14.2} {:>14.3} {f_free:>16.2}",
            v_hz * f_so * 1e3,
            v_hz * f_so / s.length()
        );
    }
    println!("   (torque per Hz {t_hz:.3e} N m; the magnet's most, m B, is {m:.3e} N m per T)");

    // -----------------------------------------------------------------------------------------
    let b5 = 5e-3;
    let f_so = m * b5 / (2.0 * PI * d);
    println!("3. Speed against the field's rate at 5 mT (step-out {f_so:.2} Hz)");
    println!(
        "   {:>8} {:>10} {:>10} {:>10}",
        "f / f_so", "field Hz", "spin Hz", "mm/s"
    );
    for k in [0.25, 0.5, 0.75, 1.0, 1.1, 1.25, 1.5, 2.0, 3.0] {
        let f = k * f_so;
        let sp = spin(f, f_so);
        println!(
            "   {k:>8.2} {f:>10.2} {sp:>10.2} {:>10.2}",
            v_hz * sp * 1e3
        );
    }

    // -----------------------------------------------------------------------------------------
    println!("4. Gravity: does it sink faster than it swims?");
    let w = (mass.mass - GLYCEROL_RHO * mass.volume) * G;
    let a: Matrix3<f64> = model.r.fixed_view::<3, 3>(0, 0).into_owned();
    let a_inv = a.try_inverse().expect("translation block invertible");
    let side = a_inv * Vector3::new(-w, 0.0, 0.0);
    let down = a_inv * Vector3::new(0.0, 0.0, -w);
    println!(
        "   net weight {:.3e} N ({:.1} mg in glycerol); centre of mass {:.2} mm ahead of the centre of buoyancy",
        w,
        w / G * 1e6,
        (mass.com - mass.cob).z * 1e3
    );
    println!(
        "   sinking, axis held level: {:.2} mm/s (drift along the axis {:+.3} mm/s)",
        -side.x * 1e3,
        side.z * 1e3
    );
    println!(
        "   sinking, axis held upright: {:.2} mm/s",
        -down.z * 1e3
    );
    let f_up = -down.z / v_hz;
    println!(
        "   to climb straight up it must turn faster than {f_up:.2} Hz; step-out at 5 mT is {f_so:.2} Hz",
    );
    // Top speed per tesla of field: v_hz * m / (2 pi D).
    let top_per_t = v_hz * m / (2.0 * PI * d);
    println!(
        "   so the field must beat {:.2} mT to outswim sinking when level, {:.2} mT to climb",
        -side.x / top_per_t * 1e3,
        -down.z / top_per_t * 1e3
    );
    let tilt = mass.mass * G * (mass.com - mass.cob).z.abs();
    println!(
        "   gravity's tilting torque, level: {tilt:.3e} N m = {:.1}% of the magnet's most at 5 mT",
        100.0 * tilt / (m * b5)
    );

    // -----------------------------------------------------------------------------------------
    println!("5. Reynolds numbers (inertia against stickiness; Stokes flow needs them far under 1)");
    let top = v_hz * f_so;
    let re = GLYCEROL_RHO * top * s.length() / GLYCEROL_MU;
    let re_rot = GLYCEROL_RHO * 2.0 * PI * f_so * h.radius.powi(2) / GLYCEROL_MU;
    println!(
        "   desk, glycerol, at step-out for 5 mT ({:.2} mm/s, {f_so:.2} Hz): Re {re:.1e} along, {re_rot:.1e} turning",
        top * 1e3
    );
    println!(
        "   rough size of the inertia error there: Oseen's first correction for a sphere as wide as the swimmer is long, (3/16) Re = {:.1}%",
        100.0 * 3.0 / 16.0 * re
    );
    for b in [1.0, 2.0, 10.0] {
        let f = m * b * 1e-3 / (2.0 * PI * d);
        let re_b = GLYCEROL_RHO * v_hz * f * s.length() / GLYCEROL_MU;
        println!("   at step-out for {b:.0} mT: Re {re_b:.1e} along");
    }
    let small = s.scaled(50e-6 / s.length());
    let f_small = 10.0;
    let u_small = advance(&h, &held) * small.helix.pitch * f_small;
    let re_s = WATER_RHO * u_small * small.length() / WATER_MU;
    let re_s_rot = WATER_RHO * 2.0 * PI * f_small * small.helix.radius.powi(2) / WATER_MU;
    println!(
        "   the same shape 50 um long in water at {f_small} Hz ({:.0} um/s): Re {re_s:.1e} along, {re_s_rot:.1e} turning",
        u_small * 1e6
    );
    let tau = mass.mass / model.r[(2, 2)];
    println!(
        "   the desk swimmer's own inertia: it settles to a new speed in {:.2} ms (mass / drag)",
        tau * 1e3
    );
}
