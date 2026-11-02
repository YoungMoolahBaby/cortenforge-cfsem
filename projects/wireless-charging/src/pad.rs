//! `pad`: a charging pad and a receiver, from the Qi reference designs (The Qi Wireless Power
//! Transfer System, Power Class 0, Part 4: Reference Designs, v1.2.3, 2017):
//!
//! - Pad: Power Transmitter design A11 (section 2.2.11): a round single-layer coil of 1.15 mm Litz,
//!   outer diameter 44 +- 1.5 mm, inner 20.5 +- 0.5 mm, 10 turns per layer; ferrite shielding at most
//!   1.0 mm under it, reaching at least 2 mm past it; the interface surface 2 mm above the coil.
//!   Spec: L_P = 6.3 uH +- 10% (coil and shielding), C_P = 0.4 uF.
//! - Receiver: Power Receiver example 1 (section 3.1): a rectangular single-layer coil of 0.41 mm
//!   Litz, outer 44.25 x 30.25 mm, inner 28.75 x 14.75 mm, 14 turns, 0.6 mm thick; ferrite
//!   52 x 35 x 1 mm on its top face; the coil 2.5 mm above the interface surface. Spec: L_s =
//!   15.3 +- 1 uH away from the pad, L'_s = 20.0 +- 1 uH on it.
//!
//! cfsem has no magnetic materials: each ferrite sheet is an ideal iron plane (images), which
//! overstates a finite sheet. With both sheets the images form a series between two planes.

use std::f64::consts::PI;

use cfsem::MU_0;
use rayon::prelude::*;

use crate::loops::{Loop, P3, mutual, self_inductance};

/// Pad coil: turn centres from r 10.825 to 21.425 mm (inner and outer diameters less the Litz
/// radius), in the plane z = 0.
pub(crate) const TX_WIRE: f64 = 0.575e-3;
const TX_TURNS: usize = 10;
/// Receiver: 14 rectangular turns at a pitch of (44.25 - 28.75) / 2 / 14 mm.
pub(crate) const RX_WIRE: f64 = 0.205e-3;
const RX_TURNS: usize = 14;
/// Planes (m): the pad's ferrite 0.15 mm (adhesive) under its coil; the receiver's coil centre
/// 2.0 (to the interface) + 2.5 + the two coils' half thicknesses above the pad coil's centre;
/// the receiver's ferrite on its top face.
pub(crate) const FERRITE_GAP: f64 = 0.15e-3;
pub(crate) const RX_Z: f64 = TX_WIRE + 2.0e-3 + 2.5e-3 + 0.3e-3;
/// Image reflections each way between two planes.
const N_IMG: usize = 6;
/// Litz resistance at ~100 kHz over DC (assumed); copper areas from the strand counts.
pub(crate) const LITZ_AC: f64 = 1.2;
pub(crate) const F: f64 = 100e3;

pub(crate) fn tx_coil(m: usize) -> Vec<Loop> {
    let (r0, r1) = (10.25e-3 + TX_WIRE, 22.0e-3 - TX_WIRE);
    (0..TX_TURNS)
        .map(|k| {
            Loop::circle(
                (0.0, 0.0, 0.0),
                r0 + (r1 - r0) * k as f64 / (TX_TURNS - 1) as f64,
                m,
                1.0,
            )
        })
        .collect()
}

pub(crate) fn rx_coil(piece: f64) -> Vec<Loop> {
    let p = (44.25e-3 - 28.75e-3) / 2.0 / RX_TURNS as f64;
    (0..RX_TURNS)
        .map(|k| {
            let s = (2 * k + 1) as f64 * p;
            Loop::rect_by(0.0, 0.0, 0.0, 44.25e-3 - s, 30.25e-3 - s, piece, 1.0)
        })
        .collect()
}

/// Images of a loop in iron planes below (a) and/or above (b).
pub(crate) fn images(l: &Loop, a: Option<f64>, b: Option<f64>) -> Vec<Loop> {
    match (a, b) {
        (None, None) => vec![],
        (Some(a), None) => vec![l.mirror(a)],
        (None, Some(b)) => vec![l.mirror(b)],
        (Some(a), Some(b)) => {
            let mut v = Vec::new();
            for first_a in [true, false] {
                let mut cur = l.clone();
                for k in 0..N_IMG {
                    let use_a = (k % 2 == 0) == first_a;
                    cur = cur.mirror(if use_a { a } else { b });
                    v.push(cur.clone());
                }
            }
            v
        }
    }
}

/// (H) series inductance of a set of turns with images.
pub(crate) fn l_coil(turns: &[Loop], wire: f64, a: Option<f64>, b: Option<f64>) -> f64 {
    turns
        .par_iter()
        .enumerate()
        .map(|(i, ti)| {
            let mut l = self_inductance(ti, wire);
            for (j, tj) in turns.iter().enumerate() {
                if i != j {
                    l += mutual(tj, ti);
                }
                for im in images(tj, a, b) {
                    l += mutual(&im, ti);
                }
            }
            l
        })
        .sum()
}

/// (H) mutual inductance, pad turns (with images) onto receiver turns.
fn m_coils(tx: &[Loop], rx: &[Loop], a: Option<f64>, b: Option<f64>) -> f64 {
    tx.par_iter()
        .map(|t| {
            let mut src = vec![t.clone()];
            src.extend(images(t, a, b));
            src.iter()
                .map(|s| rx.iter().map(|r| mutual(s, r)).sum::<f64>())
                .sum::<f64>()
        })
        .sum()
}

/// (H) Neumann's integral for two circles (centre, radius, tilt about y), trapezoid rule in both
/// angles (spectrally accurate for separate smooth loops); independent of cfsem.
fn neumann_circles(c1: P3, r1: f64, t1: f64, c2: P3, r2: f64, t2: f64, n: usize) -> f64 {
    let pt = |c: P3, r: f64, t: f64, a: f64| {
        let (x, y) = (r * a.cos(), r * a.sin());
        let (s, co) = t.sin_cos();
        (
            (co * x + c.0, y + c.1, -s * x + c.2),
            (co * -a.sin() * r, a.cos() * r, -s * -a.sin() * r),
        )
    };
    let h = 2.0 * PI / n as f64;
    let mut sum = 0.0;
    for i in 0..n {
        let (p, dp) = pt(c1, r1, t1, i as f64 * h);
        for j in 0..n {
            let (q, dq) = pt(c2, r2, t2, j as f64 * h);
            let d = ((p.0 - q.0).powi(2) + (p.1 - q.1).powi(2) + (p.2 - q.2).powi(2)).sqrt();
            sum += (dp.0 * dq.0 + dp.1 * dq.1 + dp.2 * dq.2) / d;
        }
    }
    MU_0 / (4.0 * PI) * sum * h * h
}

/// Best link efficiency at the optimal load, both sides tuned.
pub(crate) fn eta_max(k: f64, q1: f64, q2: f64) -> f64 {
    let x = k * k * q1 * q2;
    x / (1.0 + (1.0 + x).sqrt()).powi(2)
}

pub fn run() {
    println!(
        "Qi pad (A11) and receiver (example 1), Qi Part 4 v1.2.3; ferrite as ideal iron planes"
    );

    println!(
        "1. Offset and tilted circles, cfsem polygons vs an independent Neumann integral on exact circles"
    );
    for (name, c2, t2) in [
        ("offset 8 mm, 5.3 mm up", (0.008, 0.0, 0.0053), 0.0),
        (
            "offset 8 mm, tilted 10 deg",
            (0.008, 0.003, 0.0053),
            10f64.to_radians(),
        ),
    ] {
        let want = neumann_circles((0.0, 0.0, 0.0), 0.021, 0.0, c2, 0.015, t2, 400);
        for m in [64, 256, 1024] {
            let a = Loop::circle((0.0, 0.0, 0.0), 0.021, m, 1.0);
            let b = Loop::circle((0.0, 0.0, 0.0), 0.015, m, 1.0).placed(t2, c2);
            let got = mutual(&a, &b);
            println!(
                "   {name}, {m:>4}-gons: {got:.6e} vs {want:.6e} H ({:+.1e})",
                got / want - 1.0
            );
        }
    }

    // 2. Inductances against the spec.
    let a_tx = Some(-TX_WIRE - FERRITE_GAP);
    let b_rx = Some(RX_Z + 0.3e-3);
    println!(
        "2. Inductances (pad ferrite plane {:.3} mm under the pad coil, receiver ferrite on its top face)",
        (TX_WIRE + FERRITE_GAP) * 1e3
    );
    for m in [64, 128, 256] {
        let tx = tx_coil(m);
        println!(
            "   pad coil as {m}-gons: free {:.3} uH, with its ferrite {:.3} uH (10 turns in series)",
            l_coil(&tx, TX_WIRE, None, None) * 1e6,
            l_coil(&tx, TX_WIRE, a_tx, None) * 1e6
        );
    }
    let tx = tx_coil(128);
    let l_tx_free = l_coil(&tx, TX_WIRE, None, None);
    let l_tx = l_coil(&tx, TX_WIRE, a_tx, None);
    println!(
        "   spec L_P 6.3 uH +- 10%: 10 turns in series {:.2} uH ({:+.0}%); wound as 5 bifilar turns (two strands in parallel, a quarter) {:.2} uH",
        l_tx * 1e6,
        (l_tx / 6.3e-6 - 1.0) * 100.0,
        l_tx / 4.0 * 1e6
    );
    let rx_local = rx_coil(1e-3);
    let rx_away = l_coil(&rx_local, RX_WIRE, None, Some(0.3e-3));
    let rx_free = l_coil(&rx_local, RX_WIRE, None, None);
    let rx: Vec<Loop> = rx_local
        .iter()
        .map(|l| l.placed(0.0, (0.0, 0.0, RX_Z)))
        .collect();
    let rx_on = l_coil(&rx, RX_WIRE, a_tx, b_rx);
    let tx_on = l_coil(&tx, TX_WIRE, a_tx, b_rx);
    println!(
        "   receiver: free {:.2} uH; with its ferrite (away from the pad) {:.2} uH, spec 15.3 +- 1 ({:+.0}%); \
         on the pad (both sheets) {:.2} uH, spec 20.0 +- 1 ({:+.0}%)",
        rx_free * 1e6,
        rx_away * 1e6,
        (rx_away / 15.3e-6 - 1.0) * 100.0,
        rx_on * 1e6,
        (rx_on / 20.0e-6 - 1.0) * 100.0
    );
    println!(
        "   pad coil with the receiver on it (both sheets): {:.2} uH",
        tx_on * 1e6
    );

    // 3. Coupling and efficiency against offset, in free space and with both sheets.
    let len = |c: &[Loop]| c.iter().map(Loop::length).sum::<f64>();
    let r_tx = len(&tx) * 1.72e-8 / (105.0 * PI / 4.0 * 0.08e-3 * 0.08e-3) * LITZ_AC;
    let r_rx = len(&rx) * 1.72e-8 / (26.0 * PI / 4.0 * 0.08e-3 * 0.08e-3) * LITZ_AC;
    let w = 2.0 * PI * F;
    let (q1, q2) = (w * tx_on / r_tx, w * rx_on / r_rx);
    let (q1f, q2f) = (w * l_tx_free / r_tx, w * rx_free / r_rx);
    println!(
        "3. Coupling against the receiver's offset (its long side along x), {:.2} mm between the coils' centres; \
         R {:.0} / {:.0} mohm (Litz x{LITZ_AC}, assumed), Q at {:.0} kHz {:.0} / {:.0} (with ferrite), {:.0} / {:.0} (free)",
        RX_Z * 1e3,
        r_tx * 1e3,
        r_rx * 1e3,
        F / 1e3,
        q1,
        q2,
        q1f,
        q2f
    );
    // The map at 64-gons and 2 mm pieces (about 1e-3, check 1), inductances at the same.
    let (txm, rxm) = (tx_coil(64), rx_coil(2e-3));
    let rxm_on: Vec<Loop> = rxm
        .iter()
        .map(|l| l.placed(0.0, (0.0, 0.0, RX_Z)))
        .collect();
    let (lt_f, lr_f) = (
        l_coil(&txm, TX_WIRE, None, None),
        l_coil(&rxm, RX_WIRE, None, None),
    );
    let (lt_o, lr_o) = (
        l_coil(&txm, TX_WIRE, a_tx, b_rx),
        l_coil(&rxm_on, RX_WIRE, a_tx, b_rx),
    );
    let steps: Vec<f64> = (0..=10).map(|k| 2e-3 * k as f64).collect();
    let grid: Vec<(f64, f64)> = steps
        .iter()
        .flat_map(|&x| steps.iter().map(move |&y| (x, y)))
        .collect();
    let rx_at =
        |x: f64, y: f64| -> Vec<Loop> { rxm.iter().map(|l| l.placed(0.0, (x, y, RX_Z))).collect() };
    let k_map: Vec<(f64, f64)> = grid
        .par_iter()
        .map(|&(x, y)| {
            let r = rx_at(x, y);
            (
                m_coils(&txm, &r, None, None) / (lt_f * lr_f).sqrt(),
                m_coils(&txm, &r, a_tx, b_rx) / (lt_o * lr_o).sqrt(),
            )
        })
        .collect();
    let at = |x: f64, y: f64| {
        k_map[grid
            .iter()
            .position(|g| (g.0 - x).abs() < 1e-9 && (g.1 - y).abs() < 1e-9)
            .expect("grid")]
    };
    for (label, sel) in [("free space", 0usize), ("with both ferrite sheets", 1)] {
        println!("   k, {label} (rows: offset along x; columns: along y, 0 4 8 12 16 20 mm):");
        for x in [0.0, 0.004, 0.008, 0.012, 0.016, 0.02] {
            let row: Vec<String> = [0.0, 0.004, 0.008, 0.012, 0.016, 0.02]
                .iter()
                .map(|&y| {
                    let k = if sel == 0 { at(x, y).0 } else { at(x, y).1 };
                    format!("{k:+.3}")
                })
                .collect();
            println!("      x {:>2.0} mm: {}", x * 1e3, row.join("  "));
        }
        let k0 = if sel == 0 {
            at(0.0, 0.0).0
        } else {
            at(0.0, 0.0).1
        };
        let (qa, qb) = if sel == 0 { (q1f, q2f) } else { (q1, q2) };
        let cell = 4e-6 * 1e6; // mm^2 per 2 x 2 mm cell
        let area = |frac: f64| {
            // Quadrant cells (corner-weighted) where k >= frac k0, times four quadrants.
            grid.iter()
                .zip(&k_map)
                .filter(|(_, k)| (if sel == 0 { k.0 } else { k.1 }) >= frac * k0)
                .map(|(g, _)| {
                    cell * if g.0 == 0.0 { 0.5 } else { 1.0 } * if g.1 == 0.0 { 0.5 } else { 1.0 }
                })
                .sum::<f64>()
                * 4.0
        };
        println!(
            "      centred k {k0:.3}, best efficiency {:.1}%; area with k >= 90% of it {:.0} mm^2, >= 75% {:.0} mm^2, >= 50% {:.0} mm^2; \
             efficiency at 50% of k: {:.1}%",
            eta_max(k0, qa, qb) * 100.0,
            area(0.9),
            area(0.75),
            area(0.5),
            eta_max(0.5 * k0, qa, qb) * 100.0
        );
    }
}
