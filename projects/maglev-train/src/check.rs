//! `check`: forces and inductance for EMS, checked before anything is built on them.

use std::f64::consts::PI;

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_density_circular_filament;

use crate::em::{Coil, Wire, force_lorentz, image_mutual, images, lift, mutual, self_inductance};

struct Tally {
    pass: usize,
    fail: usize,
}

impl Tally {
    fn cmp(&mut self, what: &str, got: f64, want: f64, tol: f64) {
        let err = ((got - want) / want).abs();
        let ok = err <= tol;
        if ok {
            self.pass += 1;
        } else {
            self.fail += 1;
        }
        println!(
            "  {} {what:<46} {got:+.6e} vs {want:+.6e}  rel {err:.1e} (tol {tol:.0e})",
            if ok { "ok  " } else { "FAIL" }
        );
    }
}

/// K(m), E(m) by the arithmetic-geometric mean, given m and 1 - m.
fn ellip(m: f64, mc: f64) -> (f64, f64) {
    let (mut a, mut b) = (1.0_f64, mc.sqrt());
    let (mut sum, mut pow) = (0.5 * m, 0.5);
    for _ in 0..40 {
        let c = 0.5 * (a - b);
        if c.abs() < 1e-17 {
            break;
        }
        (a, b) = (0.5 * (a + b), (a * b).sqrt());
        pow *= 2.0;
        sum += pow * c * c;
    }
    let k = PI / (2.0 * a);
    (k, k * (1.0 - sum))
}

/// Maxwell's closed form for coaxial loops: M = mu0 sqrt(ab) [(2/k - k) K - (2/k) E].
fn maxwell_mutual(a: f64, b: f64, d: f64) -> f64 {
    let p2 = (a + b) * (a + b) + d * d;
    let m = 4.0 * a * b / p2;
    let k = m.sqrt();
    let (kk, ee) = ellip(m, ((a - b) * (a - b) + d * d) / p2);
    MU_0 * (a * b).sqrt() * ((2.0 / k - k) * kk - 2.0 / k * ee)
}

pub fn run() -> bool {
    let mut t = Tally { pass: 0, fail: 0 };
    let i = 10.0; // A

    println!("1. Force between coaxial loops (same current, attracting), three routes");
    println!(
        "   Lorentz (cfsem B_r) vs energy (cfsem flux, dM/dd) vs Maxwell's closed form (dM/dd)"
    );
    for (a, b, d) in [
        (0.02, 0.02, 0.001),
        (0.02, 0.02, 0.01),
        (0.02, 0.015, 0.005),
        (0.02, 0.02, 0.05),
    ] {
        let lorentz = force_lorentz(&[(a, 0.0)], i, &[(b, d)], i);
        let h = 1e-6 * d.max(1e-3);
        let m = |dd: f64| mutual(&[(a, 0.0)], &[(b, dd)]);
        let energy = -i * i * (m(d + h) - m(d - h)) / (2.0 * h);
        let closed =
            -i * i * (maxwell_mutual(a, b, d + h) - maxwell_mutual(a, b, d - h)) / (2.0 * h);
        t.cmp(
            &format!("a {a}, b {b}, d {d} m: Lorentz vs closed form"),
            lorentz,
            closed,
            1e-6,
        );
        t.cmp("    energy vs closed form", energy, closed, 1e-6);
    }
    // Close loops look like parallel wires: F = mu0 I^2 a / d, corrections O((d/a) ln(a/d)).
    let (a, d) = (0.02, 2e-5);
    t.cmp(
        "d/a = 1e-3: Lorentz vs parallel wires mu0 I^2 a/d",
        force_lorentz(&[(a, 0.0)], i, &[(a, d)], i),
        MU_0 * i * i * a / d,
        2e-2,
    );

    println!(
        "2. The image across an iron surface: tangential B vanishes on the plane (same-sign image)"
    );
    let lp = [(0.02, -0.005)];
    let src = [lp[0], images(&lp)[0]];
    let (rs, zs): (Vec<f64>, Vec<f64>) = src.iter().copied().unzip();
    for r in [0.01, 0.02, 0.04] {
        let (mut br, mut bz) = ([0.0], [0.0]);
        flux_density_circular_filament(
            (&rs, &zs, &[1.0, 1.0]),
            &[0.0, 0.0],
            (&[r], &[0.0]),
            (&mut br, &mut bz),
        )
        .expect("lengths");
        let ok = br[0].abs() <= 1e-12 * bz[0].abs();
        if ok {
            t.pass += 1
        } else {
            t.fail += 1
        }
        println!(
            "  {} at r {r} m on the plane: B_r {:+.2e} T, B_z {:+.4e} T (B_r / B_z {:.1e})",
            if ok { "ok  " } else { "FAIL" },
            br[0],
            bz[0],
            (br[0] / bz[0]).abs()
        );
    }

    println!("3. A pancake coil (AWG 24, r_in 6 mm, 24 layers x 17 turns) under the rail");
    let coil = Coil {
        r_in: 0.006,
        layers: 24,
        per_layer: 17,
        wire: Wire::awg(24),
    };
    println!(
        "   {} turns of AWG {}, r_out {:.2} mm, height {:.2} mm, R {:.3} ohm",
        coil.turns(),
        coil.wire.awg,
        coil.r_out() * 1e3,
        coil.height() * 1e3,
        coil.resistance()
    );
    for g in [0.002, 0.005, 0.01] {
        let f = lift(&coil, g, 1.0);
        let h = 1e-6;
        let energy = -0.5 * (image_mutual(&coil, g + h) - image_mutual(&coil, g - h)) / (2.0 * h);
        t.cmp(
            &format!("lift at gap {} mm, 1 A: Lorentz vs coenergy", g * 1e3),
            f,
            energy,
            1e-6,
        );
    }
    let l0 = self_inductance(&coil);
    // Wheeler (1928) multilayer: L[uH] = 0.8 a^2 N^2 / (6a + 9b + 10c), inches.
    let inch = 0.0254;
    let (am, b, c) = (
        0.5 * (coil.r_in + coil.r_out()) / inch,
        coil.height() / inch,
        (coil.r_out() - coil.r_in) / inch,
    );
    let n = coil.turns() as f64;
    let wheeler = 0.8 * am * am * n * n / (6.0 * am + 9.0 * b + 10.0 * c) * 1e-6;
    t.cmp(
        "free-space inductance vs Wheeler (quoted ~1%; c/a large here)",
        l0,
        wheeler,
        0.05,
    );

    println!(
        "4. A disc magnet and its image, as current loops, vs two point dipoles (cfsem's dipole field)"
    );
    {
        // 12 x 4 mm N42 at a 6 mm gap; far away the loops must look like dipoles of m = Br V / mu0.
        let (rad, h, br, g) = (0.006, 0.004, 1.3, 0.006);
        let (m, im) = crate::em::magnet_loops(rad, h, br, 80, g);
        let src: Vec<(f64, f64)> = m.iter().copied().chain(images(&m)).collect();
        let (rs, zs): (Vec<f64>, Vec<f64>) = src.iter().copied().unzip();
        let moment = br / MU_0 * PI * rad * rad * h;
        let zc = -g - 0.5 * h;
        for (x, z) in [(0.1, -0.1), (0.3, -0.2), (0.0, -1.0)] {
            let (mut brr, mut bzz) = ([0.0], [0.0]);
            flux_density_circular_filament(
                (&rs, &zs, &vec![im; src.len()]),
                &vec![0.0; src.len()],
                (&[x], &[z]),
                (&mut brr, &mut bzz),
            )
            .expect("lengths");
            let loops = brr[0].hypot(bzz[0]);
            let mut d = (0.0, 0.0, 0.0);
            for zz in [zc, -zc] {
                let b = cfsem::physics::point_source::dipole::flux_density_dipole_scalar(
                    (0.0, 0.0, zz),
                    (0.0, 0.0, moment),
                    0.0,
                    (x, 0.0, z),
                );
                d = (d.0 + b.0, d.1 + b.1, d.2 + b.2);
            }
            let dipole = d.0.hypot(d.1).hypot(d.2);
            // Corrections are O((size / distance)^2); size ~ the image pair's extent, 2 (g + h).
            let tol = 3.0 * (2.0 * (g + h) / x.hypot(z)).powi(2);
            t.cmp(
                &format!("|B| at ({x}, {z}) m: loops vs dipoles"),
                loops,
                dipole,
                tol,
            );
        }
    }

    println!(
        "5. Design numbers (not tests): what air-core EMS costs, and a permanent-magnet hybrid"
    );
    let r = coil.resistance();
    let f1 = lift(&coil, 0.005, 1.0);
    let need = 0.3 * 9.81 / 4.0;
    let i_need = (need / f1).sqrt();
    println!(
        "   air-core coil at 5 mm: {f1:.4} N at 1 A; a 300 g vehicle on 4 coils needs {i_need:.2} A each, \
         {:.0} W in total",
        4.0 * i_need * i_need * r
    );
    let (rad, h, br) = (0.010, 0.005, 1.3); // 20 x 5 mm N42 disc
    // Convergence of the sheet-current model in the number of loops.
    for n in [10, 40, 160] {
        let (m, im) = crate::em::magnet_loops(rad, h, br, n, 0.005);
        println!(
            "   N42 disc 20 x 5 mm (Br 1.3 T) at 5 mm, {n:>3} loops: lift {:.4} N",
            force_lorentz(&m, im, &images(&m), im)
        );
    }
    println!(
        "   disc lift vs gap, and the field on the rail surface (saturation check, steel ~1.5 T):"
    );
    for g in [0.002, 0.005, 0.008, 0.012] {
        let (m, im) = crate::em::magnet_loops(rad, h, br, 80, g);
        let f = force_lorentz(&m, im, &images(&m), im);
        let src: Vec<(f64, f64)> = m.iter().copied().chain(images(&m)).collect();
        let (rs, zs): (Vec<f64>, Vec<f64>) = src.iter().copied().unzip();
        let pts: Vec<f64> = (0..40).map(|k| 0.0005 * k as f64).collect();
        let (mut brr, mut bzz) = (vec![0.0; pts.len()], vec![0.0; pts.len()]);
        flux_density_circular_filament(
            (&rs, &zs, &vec![im; src.len()]),
            &vec![0.0; src.len()],
            (&pts, &vec![0.0; pts.len()]),
            (&mut brr, &mut bzz),
        )
        .expect("lengths");
        let bmax = bzz.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        println!(
            "     gap {:>4.1} mm: lift {f:.3} N ({:.0} g), max |B| on the rail {bmax:.3} T",
            g * 1e3,
            f / 9.81 * 1e3
        );
    }
    // Control authority: a coil around the disc (same axis, top flush with the magnet) adds a
    // force linear in its current, dF/dI = cross term of (magnet + coil) with their images.
    let g = 0.005;
    let ctrl = Coil {
        r_in: 0.0105,
        layers: 8,
        per_layer: 9,
        wire: Wire::awg(24),
    };
    let (m, im) = crate::em::magnet_loops(rad, h, br, 80, g);
    let c = ctrl.turns_at(g);
    let total = |ic: f64| {
        // All loops with their currents, against all images.
        let all: Vec<((f64, f64), f64)> = m
            .iter()
            .map(|&p| (p, im))
            .chain(c.iter().map(|&p| (p, ic)))
            .collect();
        let mut f = 0.0;
        for &(p, ip) in &all {
            for &(q, iq) in &all {
                f += force_lorentz(&[p], ip, &images(&[q]), iq);
            }
        }
        f
    };
    let (f0, fp) = (total(0.0), total(0.5));
    println!(
        "   control coil around the disc ({} turns AWG 24, {:.2} ohm): {:+.4} N per amp at 5 mm \
         (pull {f0:.3} N at 0 A, {fp:.3} N at +0.5 A)",
        ctrl.turns(),
        ctrl.resistance(),
        (fp - f0) / 0.5
    );

    println!("\n{} passed, {} failed", t.pass, t.fail);
    t.fail == 0
}
