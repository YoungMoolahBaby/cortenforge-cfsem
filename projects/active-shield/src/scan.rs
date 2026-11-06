//! `transmission`: the fraction of cosmic rays of each rigidity that reach the habitat through
//! the shield, relative to the habitat's view with the field off (the ring shades part of its
//! own sky), overall and at the outer and inner equators; and how that cutoff moves with the
//! bending power and the shell's thickness.

use std::time::Instant;

use crate::design;
use crate::shield::Shield;
use crate::trace::{transmission as tx, Tracer, C};

pub const RIGIDITIES_GV: [f64; 16] = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 200.0];
pub const POLOIDAL: usize = 24;
pub const DIRECTIONS: usize = 400;

/// Transmission against rigidity, normalised by the field-off view. Returns (GV, overall,
/// outer equator, inner equator).
pub fn curve(sh: &Shield, tr: &Tracer, gvs: &[f64]) -> Vec<(f64, f64, f64, f64)> {
    let off = Shield { amp_turns: 0.0, ..sh.clone() };
    let (open, open_rows) = tx(&off, tr, 1e9, POLOIDAL, DIRECTIONS);
    gvs.iter()
        .map(|&gv| {
            let (total, rows) = tx(sh, tr, gv * 1e9, POLOIDAL, DIRECTIONS);
            let outer = rows[0].1 / open_rows[0].1;
            let inner = rows[POLOIDAL / 2].1 / open_rows[POLOIDAL / 2].1;
            (gv, total / open, outer, inner)
        })
        .collect()
}

/// Rigidity (GV) where the transmission first reaches `level`, by linear interpolation.
pub fn rigidity_at(curve: &[(f64, f64, f64, f64)], level: f64) -> f64 {
    for w in curve.windows(2) {
        if w[0].1 < level && w[1].1 >= level {
            return w[0].0 + (level - w[0].1) / (w[1].1 - w[0].1) * (w[1].0 - w[0].0);
        }
    }
    f64::NAN
}

pub fn transmission() {
    let sh = design::baseline();
    let tr = Tracer { r_hab: design::R_HAB, ..Tracer::default() };
    let off = Shield { amp_turns: 0.0, ..sh.clone() };
    let (open, rows) = tx(&off, &tr, 1e9, POLOIDAL, DIRECTIONS);
    println!(
        "Baseline: {:.0} T m ({:.1} T at R = {:.0} m, shell {:.0}-{:.0} m), {} poloidal points x {} directions",
        sh.bending_power(), sh.b_at(sh.r_major), sh.r_major, sh.r1, sh.r2, POLOIDAL, DIRECTIONS
    );
    println!(
        "  field off: the habitat sees {:.4} of the sky (outer equator {:.4}, inner {:.4}; the ring shades the rest)",
        open, rows[0].1, rows[POLOIDAL / 2].1
    );
    println!("  c B d at the tube centre's radius: {:.1} GV\n", C * sh.bending_power() / 1e9);
    println!("  {:>6}  {:>8}  {:>8}  {:>8}", "GV", "overall", "outer eq", "inner eq");
    let t0 = Instant::now();
    let base = curve(&sh, &tr, &RIGIDITIES_GV);
    for (gv, t, o, i) in &base {
        println!("  {gv:>6.0}  {t:>8.4}  {o:>8.4}  {i:>8.4}");
    }
    println!("  ({:.1} s)", t0.elapsed().as_secs_f64());
    println!(
        "  half-transmission rigidity {:.1} GV, 10% at {:.1} GV, 90% at {:.1} GV\n",
        rigidity_at(&base, 0.5),
        rigidity_at(&base, 0.1),
        rigidity_at(&base, 0.9)
    );

    println!("Bending power and thickness (inner winding at {:.0} m): rigidity (GV) where 10% / 50% / 90% get through", sh.r1);
    println!("  {:>6}  {:>5}  {:>6}  {:>7}  {:>6}  {:>6}  {:>6}  {:>10}", "T m", "d m", "B T", "c B d", "10%", "50%", "90%", "energy J");
    for (bl, d) in [(25.0, 20.0), (50.0, 20.0), (100.0, 10.0), (100.0, 20.0), (100.0, 40.0), (200.0, 20.0), (200.0, 40.0)] {
        let b = bl / d;
        let s = Shield::with_field(sh.r_major, sh.r1, sh.r1 + d, sh.n, b);
        let c = curve(&s, &tr, &RIGIDITIES_GV);
        println!(
            "  {bl:>6.0}  {d:>5.0}  {b:>6.2}  {:>7.1}  {:>6.1}  {:>6.1}  {:>6.1}  {:>10.3e}",
            C * bl / 1e9,
            rigidity_at(&c, 0.1),
            rigidity_at(&c, 0.5),
            rigidity_at(&c, 0.9),
            s.energy_ideal()
        );
    }
}
