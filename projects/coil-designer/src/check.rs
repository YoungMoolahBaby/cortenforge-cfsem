//! `check`: cfsem's circular-filament API against known results, before anything is built on it.
//!
//! Each line compares a cfsem value with a reference that does not come from cfsem: a textbook
//! closed form, a series, or a quadrature written here.

use std::f64::consts::PI;
use std::time::Instant;

use cfsem::MU_0;

use crate::coil::{Coil, Pack, Wire, ring_self_inductance, ring_self_inductance_n, uniform_radius};
use crate::refs::{
    biot_savart_loop, closed_form_loop, current_sheet, kirchhoff_ring, maxwell_mutual,
    neumann_coaxial, thick_solenoid_centre, wheeler_multilayer, wheeler_single,
};
use cfsem::physics::circular_filament::{
    flux_circular_filament, flux_density_circular_filament, flux_density_circular_filament_on_axis,
    flux_density_circular_filament_par, flux_density_circular_filament_scalar,
};

/// Running tally of checks.
struct Tally {
    pass: usize,
    fail: usize,
}

impl Tally {
    /// Print one comparison and count it. `tol` is on the relative error.
    fn cmp(&mut self, what: &str, got: f64, want: f64, tol: f64) {
        let err = rel(got, want);
        self.verdict(
            what,
            &format!("{got:+.9e} vs {want:+.9e}  rel {err:.1e}"),
            err <= tol,
            tol,
        );
    }

    /// Print a precomputed error and count it.
    fn err(&mut self, what: &str, detail: &str, err: f64, tol: f64) {
        self.verdict(what, &format!("{detail}  rel {err:.1e}"), err <= tol, tol);
    }

    fn verdict(&mut self, what: &str, detail: &str, ok: bool, tol: f64) {
        if ok {
            self.pass += 1;
        } else {
            self.fail += 1;
        }
        let tag = if ok { "ok  " } else { "FAIL" };
        println!("  {tag} {what:<34} {detail}  (tol {tol:.0e})");
    }
}

fn rel(got: f64, want: f64) -> f64 {
    if want == 0.0 {
        got.abs()
    } else {
        (got - want).abs() / want.abs()
    }
}

/// Sum (Br, Bz) of several ideal loops `(r, z, NI)` at one point, through cfsem's vector API.
fn field(loops: &[(f64, f64, f64)], (r, z): (f64, f64)) -> (f64, f64) {
    let rf: Vec<f64> = loops.iter().map(|l| l.0).collect();
    let zf: Vec<f64> = loops.iter().map(|l| l.1).collect();
    let nif: Vec<f64> = loops.iter().map(|l| l.2).collect();
    let wire = vec![0.0; loops.len()];
    let (mut br, mut bz) = ([0.0], [0.0]);
    flux_density_circular_filament((&rf, &zf, &nif), &wire, (&[r], &[z]), (&mut br, &mut bz))
        .expect("lengths match");
    (br[0], bz[0])
}

/// Checks 9-12: the turn-by-turn coil model in `coil.rs`.
fn winding_packs(t: &mut Tally) {
    let wire = Wire::awg(22);
    println!("9. Winding-pack inductance: cfsem kernels vs Kirchhoff + Maxwell closed form");
    let pack = Pack {
        r_in: 0.02,
        z_mid: 0.0,
        layers: 4,
        per_layer: 6,
        wire,
        sign: 1.0,
    };
    let coil = Coil { packs: vec![pack] };
    let f = coil.filaments();
    let mut want = 0.0;
    for (i, &(ri, zi, _, bi)) in f.iter().enumerate() {
        want += kirchhoff_ring(ri, bi);
        for &(rj, zj, _, _) in &f[i + 1..] {
            want += 2.0 * maxwell_mutual(ri, rj, zj - zi);
        }
    }
    // Kirchhoff drops O((b/a)^2) terms; here b/a is about 0.015.
    t.cmp(
        "4 x 6 turns, AWG 22, r_in 20 mm",
        coil.inductance(),
        want,
        1e-4,
    );

    println!("10. Multilayer coil vs Wheeler (1928), quoted to ~1%");
    // Wheeler is most accurate when 6a, 9b and 10c are comparable: a = 1 in, b ~ 2/3 in, c ~ 0.6 in.
    let (layers, per_layer) = (22, 24);
    let depth = layers as f64 * wire.outer;
    let pack = Pack {
        r_in: 0.0254 - 0.5 * depth,
        z_mid: 0.0,
        layers,
        per_layer,
        wire,
        sign: 1.0,
    };
    let coil = Coil { packs: vec![pack] };
    let want = wheeler_multilayer(
        pack.r_mean(),
        pack.width(),
        pack.depth(),
        pack.turns() as f64,
    );
    t.cmp(
        &format!("{} turns, 1 in mean radius", pack.turns()),
        coil.inductance(),
        want,
        0.02,
    );

    println!(
        "11. Single-layer solenoid: current sheet (Nagaoka, numeric) vs Wheeler and filaments"
    );
    let per_layer = 72;
    let pack = Pack {
        r_in: 0.0254 - 0.5 * wire.outer,
        z_mid: 0.0,
        layers: 1,
        per_layer,
        wire,
        sign: 1.0,
    };
    let (a, l, n) = (pack.r_mean(), pack.width(), per_layer as f64);
    let sheet = current_sheet(a, l, n);
    t.cmp(
        "Wheeler vs sheet (l = 2a)",
        wheeler_single(a, l, n),
        sheet,
        0.01,
    );
    let fil = Coil { packs: vec![pack] }.inductance();
    // Round wire stores less energy than a uniform sheet (Rosa's correction): expect a small
    // negative difference, of order mu0 a N / L_sheet.
    let scale = MU_0 * a * n / sheet;
    println!(
        "  info filaments {fil:.4e} H, sheet {sheet:.4e} H, diff {:+.2e} (mu0 a N / L = {scale:.2e})",
        fil / sheet - 1.0
    );

    println!("12. Thick-solenoid centre field vs Montgomery's uniform-current-density formula");
    for (layers, per_layer) in [(10, 10), (30, 8), (5, 60)] {
        let pack = Pack {
            r_in: 0.03,
            z_mid: 0.0,
            layers,
            per_layer,
            wire,
            sign: 1.0,
        };
        let coil = Coil { packs: vec![pack] };
        let got = coil.field(1.0, &[0.0], &[0.0]).1[0];
        let want = thick_solenoid_centre(
            pack.r_in,
            pack.r_in + pack.depth(),
            pack.width(),
            pack.turns() as f64,
        );
        // Turns at cell centres are a midpoint rule over the section: O(pitch^2) error.
        t.cmp(&format!("{layers} layers x {per_layer}"), got, want, 1e-3);
    }

    println!("12b. Uniform sphere of a thin Helmholtz pair vs the leading-order series");
    // B/B0 - 1 = -(144/125)(z/R)^4 on the axis, the worst direction (P4 peaks there), so
    // the eps sphere has radius (125 eps / 144)^(1/4) R plus O((rho/R)^2) corrections.
    let a = 0.05;
    let helm = [(a, -a / 2.0, 1.0), (a, a / 2.0, 1.0)];
    let f = |r: &[f64], z: &[f64]| {
        let (mut br, mut bz) = (Vec::new(), Vec::new());
        for (&r, &z) in r.iter().zip(z) {
            let (x, y) = field(&helm, (r, z));
            br.push(x);
            bz.push(y);
        }
        (br, bz)
    };
    for eps in [1e-6, 1e-8] {
        let want = (125.0 * eps / 144.0_f64).powf(0.25) * a;
        t.cmp(
            &format!("eps {eps:.0e}"),
            uniform_radius(&f, eps, a),
            want,
            1e-2,
        );
    }
}

pub fn run() -> bool {
    let mut t = Tally { pass: 0, fail: 0 };
    let a = 0.1; // m, loop radius
    let ni = 100.0; // A-turns

    // cfsem's elliptic integrals are curve fits (math.rs: Abramowitz & Stegun 17.3.34/36, max
    // absolute error 2e-8; its D fit, 4e-9 relative), so its fields carry ~1e-8 relative error.
    const FIT: f64 = 3e-8;

    println!("1. Single loop on the axis: Bz = mu0 I a^2 / (2 (a^2 + z^2)^1.5)");
    for zr in [0.0, 0.5, 2.0, 10.0] {
        let z = zr * a;
        let want = MU_0 * ni * a * a / (2.0 * (a * a + z * z).powf(1.5));
        let got = flux_density_circular_filament_scalar((a, 0.0, ni), (0.0, z)).1;
        t.cmp(&format!("z = {zr} a"), got, want, FIT);
        let got = flux_density_circular_filament_on_axis((a, 0.0, ni), z).1;
        t.cmp(&format!("z = {zr} a (_on_axis)"), got, want, 1e-14);
    }

    println!("2. Off-axis field: cfsem and an exact closed form (AGM) vs Biot-Savart quadrature");
    for (rr, zr) in [
        (0.5, 0.0),
        (0.5, 0.3),
        (0.9, 0.05),
        (1.05, -0.02),
        (2.0, 1.0),
        (3.0, -2.0),
    ] {
        let p = (rr * a, zr * a);
        let (qr, qz) = biot_savart_loop(a, ni, p, 1 << 16);
        let mag = qr.hypot(qz);
        let (br, bz) = flux_density_circular_filament_scalar((a, 0.0, ni), p);
        let err = (br - qr).hypot(bz - qz) / mag;
        t.err(
            &format!("cfsem at ({rr}, {zr}) a"),
            &format!("|B| {mag:.6e} T"),
            err,
            FIT,
        );
        // The reference itself: quadrature and exact closed form must agree to rounding.
        let (cr, cz) = closed_form_loop(a, ni, p);
        let err = (cr - qr).hypot(cz - qz) / mag;
        t.err(
            &format!("  AGM closed form at ({rr}, {zr}) a"),
            "",
            err,
            1e-13,
        );
    }

    println!("3. Far field against a point dipole, m = I pi a^2");
    let m = ni * PI * a * a;
    for dr in [100.0, 1e4] {
        let d = dr * a;
        let k = MU_0 * m / (4.0 * PI * d.powi(3));
        let ax = flux_density_circular_filament_scalar((a, 0.0, ni), (0.0, d)).1;
        // The dipole limit has O((a/d)^2) corrections.
        let tol = 3.0 / (dr * dr);
        t.cmp(&format!("axis, z = {dr:.0e} a"), ax, 2.0 * k, tol);
        let eq = flux_density_circular_filament_scalar((a, 0.0, ni), (d, 0.0)).1;
        t.cmp(&format!("equator, r = {dr:.0e} a"), eq, -k, tol);
    }

    println!("4. Maxwell's equations by central differences (h = 1e-5 a), ideal loop");
    for (rr, zr) in [(0.3, 0.2), (1.2, 0.1), (2.0, -1.5)] {
        let (r, z, h) = (rr * a, zr * a, 1e-5 * a);
        let b = |r, z| flux_density_circular_filament_scalar((a, 0.0, ni), (r, z));
        // div B = (1/r) d(r Br)/dr + dBz/dz; curl_phi B = dBr/dz - dBz/dr.
        let drbr = ((r + h) * b(r + h, z).0 - (r - h) * b(r - h, z).0) / (2.0 * h);
        let dbzdz = (b(r, z + h).1 - b(r, z - h).1) / (2.0 * h);
        let dbrdz = (b(r, z + h).0 - b(r, z - h).0) / (2.0 * h);
        let dbzdr = (b(r + h, z).1 - b(r - h, z).1) / (2.0 * h);
        let scale = b(r, z).0.hypot(b(r, z).1) / a; // T/m, the size of a typical derivative
        let div = (drbr / r + dbzdz).abs() / scale;
        let curl = (dbrdz - dbzdr).abs() / scale;
        t.err(&format!("div B at ({rr}, {zr}) a"), "", div, 1e-6);
        t.err(&format!("curl B at ({rr}, {zr}) a"), "", curl, 1e-6);
    }

    println!("5. Helmholtz pair (spacing = R): centre field and 4th-order flatness");
    let helm = [(a, -a / 2.0, ni), (a, a / 2.0, ni)];
    let b0 = field(&helm, (0.0, 0.0)).1;
    t.cmp(
        "Bz(0) = (4/5)^1.5 mu0 NI / R",
        b0,
        (0.8_f64).powf(1.5) * MU_0 * ni / a,
        FIT,
    );
    for zr in [0.02, 0.05] {
        // On the axis B(z)/B0 = 1 - (144/125)(z/R)^4 + O((z/R)^6).
        let dev = field(&helm, (0.0, zr * a)).1 / b0 - 1.0;
        let series = -144.0 / 125.0 * zr.powi(4);
        // The next term is O(zr^2) relative to this one.
        t.cmp(
            &format!("B(z)/B0 - 1 at z = {zr} R"),
            dev,
            series,
            2.0 * zr * zr,
        );
    }

    println!("6. Maxwell gradient pair (spacing = sqrt(3) R, opposed currents)");
    let d = 3f64.sqrt() * a;
    let maxwell = [(a, -d / 2.0, -ni), (a, d / 2.0, ni)];
    let h = 1e-5 * a;
    let g = (field(&maxwell, (0.0, h)).1 - field(&maxwell, (0.0, -h)).1) / (2.0 * h);
    let g_text = 3.0 * MU_0 * ni * a * a * (d / 2.0) / (a * a + d * d / 4.0).powf(2.5);
    t.cmp("dBz/dz(0) closed form", g, g_text, 1e-8);
    t.cmp(
        "dBz/dz(0) = 0.6413 mu0 NI / R^2",
        g,
        0.6413 * MU_0 * ni / (a * a),
        1e-4,
    );
    // At sqrt(3) R the cubic term vanishes, so Bz(z) / (G z) - 1 = O((z/R)^4).
    let lin = field(&maxwell, (0.0, 0.1 * a)).1 / (g * 0.1 * a) - 1.0;
    t.err(
        "linearity at z = 0.1 R",
        &format!("Bz/(G z) - 1 = {lin:+.2e}"),
        lin.abs(),
        1e-3,
    );

    println!("7. Mutual inductance of coaxial loops: cfsem flux vs Neumann quadrature");
    for (b, gap) in [(0.1, 0.05), (0.05, 0.0), (0.2, 0.3), (0.1, 0.002)] {
        let mut psi = [0.0];
        flux_circular_filament((&[a], &[0.0], &[1.0]), &[0.0], (&[b], &[gap]), &mut psi)
            .expect("lengths match");
        let want = neumann_coaxial(a, b, gap, 1 << 16);
        t.cmp(&format!("a = 0.1, b = {b}, d = {gap} m"), psi[0], want, FIT);
        // The closed form used by the winding checks below, against the same quadrature.
        t.cmp(
            "  Maxwell closed form (AGM)",
            maxwell_mutual(a, b, gap),
            want,
            1e-12,
        );
    }

    println!("8. Thin-ring self-inductance: cross-section mean of cfsem's finite-wire flux");
    println!("   vs Kirchhoff, L = mu0 a (ln(8a/b) - 7/4) for uniform current, b << a");
    for ba in [0.03, 0.01, 0.001] {
        let b = ba * a;
        // L = (1/I) <flux> over the conductor section (energy method, uniform J).
        let l = ring_self_inductance(a, b);
        // The formula drops O((b/a)^2 ln) terms.
        t.cmp(&format!("b/a = {ba}"), l, kirchhoff_ring(a, b), 1e-3);
        // The quadrature is exact for cfsem's interior model: more angular points change nothing.
        t.cmp(
            "  16 vs 256 angular points",
            l,
            ring_self_inductance_n(a, b, 256),
            1e-13,
        );
    }

    winding_packs(&mut t);

    println!("13. Speed: Helmholtz field on a 1000 x 1000 (r, z) grid");
    let n = 1000;
    let (mut rs, mut zs) = (Vec::with_capacity(n * n), Vec::with_capacity(n * n));
    for i in 0..n {
        for j in 0..n {
            rs.push(2.0 * a * (i as f64 + 0.5) / n as f64);
            zs.push(-1.5 * a + 3.0 * a * (j as f64 + 0.5) / n as f64);
        }
    }
    let (mut br, mut bz) = (vec![0.0; n * n], vec![0.0; n * n]);
    let src = (&[a, a][..], &[-a / 2.0, a / 2.0][..], &[ni, ni][..]);
    let t0 = Instant::now();
    flux_density_circular_filament(src, &[0.0, 0.0], (&rs, &zs), (&mut br, &mut bz))
        .expect("lengths match");
    let serial = t0.elapsed();
    let t0 = Instant::now();
    flux_density_circular_filament_par(src, &[0.0, 0.0], (&rs, &zs), (&mut br, &mut bz))
        .expect("lengths match");
    let par = t0.elapsed();
    let per = serial.as_secs_f64() * 1e9 / (2 * n * n) as f64;
    println!(
        "  {} points x 2 loops: serial {:.1} ms ({per:.1} ns per loop-point), parallel {:.1} ms",
        n * n,
        serial.as_secs_f64() * 1e3,
        par.as_secs_f64() * 1e3
    );

    println!("\n{} passed, {} failed", t.pass, t.fail);
    t.fail == 0
}
