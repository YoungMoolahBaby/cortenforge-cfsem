//! `pack`: a Helmholtz pair wound from real wire, against the ideal thin-loop pair.

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_density_circular_filament_par;

use crate::coil::{Coil, Pack, Wire, uniform_radius};

/// Field of ideal loops `(r, z, NI)` at points, through cfsem.
pub fn thin_field(loops: &[(f64, f64, f64)], r: &[f64], z: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let rf: Vec<f64> = loops.iter().map(|l| l.0).collect();
    let zf: Vec<f64> = loops.iter().map(|l| l.1).collect();
    let nif: Vec<f64> = loops.iter().map(|l| l.2).collect();
    let (mut br, mut bz) = (vec![0.0; r.len()], vec![0.0; r.len()]);
    flux_density_circular_filament_par(
        (&rf, &zf, &nif),
        &vec![0.0; loops.len()],
        (r, z),
        (&mut br, &mut bz),
    )
    .expect("lengths match");
    (br, bz)
}

/// Two identical packs, mean radius `r_mean`, centres `spacing` apart, same current sense.
pub fn helmholtz(r_mean: f64, spacing: f64, layers: usize, per_layer: usize, wire: Wire) -> Coil {
    let r_in = r_mean - 0.5 * layers as f64 * wire.outer;
    let pack = |z_mid| Pack {
        r_in,
        z_mid,
        layers,
        per_layer,
        wire,
        sign: 1.0,
        partial: 0,
    };
    Coil {
        packs: vec![pack(-0.5 * spacing), pack(0.5 * spacing)],
    }
}

pub fn run() {
    let (r_mean, layers, per_layer, wire) = (0.05, 10, 10, Wire::awg(22));
    let coil = helmholtz(r_mean, r_mean, layers, per_layer, wire);
    let p = coil.packs[0];
    let n = p.turns() as f64;

    println!("Helmholtz pair from AWG {} wire", wire.awg);
    println!(
        "  wire: {:.3} mm copper, {:.3} mm over enamel (assumed {:.2} mm build), {:.1} mohm/m",
        wire.bare * 1e3,
        wire.outer * 1e3,
        (wire.outer - wire.bare) * 1e3,
        wire.ohm_per_m() * 1e3
    );
    println!(
        "  each coil: {} turns ({} layers x {}), radius {:.2}-{:.2} mm (mean {:.1}), width {:.2} mm",
        p.turns(),
        layers,
        per_layer,
        p.r_in * 1e3,
        (p.r_in + p.depth()) * 1e3,
        p.r_mean() * 1e3,
        p.width() * 1e3
    );
    println!(
        "  spacing: {:.1} mm between pack centres (= mean radius)",
        r_mean * 1e3
    );

    let b_per_a = coil.field(1.0, &[0.0], &[0.0]).1[0];
    let r = coil.resistance();
    let l = coil.inductance();
    println!("Electrical (series, 20 C)");
    println!(
        "  wire {:.1} m, copper {:.0} g",
        coil.wire_length(),
        coil.copper_mass() * 1e3
    );
    println!(
        "  R = {r:.3} ohm, L = {:.3} mH, L/R = {:.3} ms",
        l * 1e3,
        l / r * 1e3
    );
    let thin_b = (0.8_f64).powf(1.5) * MU_0 * n / r_mean;
    println!(
        "  B0 = {:.4} mT/A (thin-loop formula at the mean radius: {:.4} mT/A, {:+.2e})",
        b_per_a * 1e3,
        thin_b * 1e3,
        b_per_a / thin_b - 1.0
    );
    let i1 = 1e-3 / b_per_a;
    println!(
        "  1 mT takes {i1:.3} A, {:.2} V, {:.2} W",
        i1 * r,
        i1 * i1 * r
    );

    println!("Uniform sphere radius, |B - B0| / B0 <= eps (fraction of the mean radius)");
    let thin = [(r_mean, -0.5 * r_mean, n), (r_mean, 0.5 * r_mean, n)];
    let thin_f = |r: &[f64], z: &[f64]| thin_field(&thin, r, z);
    let real_f = |r: &[f64], z: &[f64]| coil.field(1.0, r, z);
    for eps in [1e-2, 1e-3, 1e-4, 1e-5] {
        let t = uniform_radius(&thin_f, eps, r_mean);
        let w = uniform_radius(&real_f, eps, r_mean);
        println!(
            "  eps {eps:.0e}: thin loops {:.2} mm ({:.3} R), real packs {:.2} mm ({:.3} R)",
            t * 1e3,
            t / r_mean,
            w * 1e3,
            w / r_mean
        );
    }

    println!("Spacing scan for the real packs (eps 1e-4)");
    let mut best = (0.0, 0.0);
    for k in 0..=40 {
        let d = r_mean * (0.95 + 0.0025 * k as f64);
        let c = helmholtz(r_mean, d, layers, per_layer, wire);
        let f = |r: &[f64], z: &[f64]| c.field(1.0, r, z);
        let u = uniform_radius(&f, 1e-4, r_mean);
        if u > best.1 {
            best = (d, u);
        }
    }
    println!(
        "  best spacing {:.4} R ({:.2} mm): 1e-4 sphere {:.2} mm ({:.3} R)",
        best.0 / r_mean,
        best.0 * 1e3,
        best.1 * 1e3,
        best.1 / r_mean
    );
}
