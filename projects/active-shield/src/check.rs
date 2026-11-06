//! `check`: cfsem's field of the discrete windings, their stored energy and forces, and the
//! particle tracer, against closed forms before any design is built on them.

use std::f64::consts::PI;

use cfsem::MU_0;

use crate::design;
use crate::dose;
use crate::shield::{cross, dot, loop_flux, loop_point, norm, Loops, Shield};
use crate::trace::{boris_turn, brho, Fate, Tracer, C};

struct Tally {
    pass: usize,
    fail: usize,
}

impl Tally {
    fn err(&mut self, what: &str, detail: &str, err: f64, tol: f64) {
        let ok = err <= tol;
        if ok {
            self.pass += 1;
        } else {
            self.fail += 1;
        }
        let tag = if ok { "ok  " } else { "FAIL" };
        println!("  {tag} {what:<48} {detail} err {err:.2e} (tol {tol:.0e})");
    }
}

/// Total radial (minor-radius) force on the loop of radius r at phi = 0, the integral of
/// (I dl x B) . rho_hat, and the force vector.
fn loop_forces(loops: &Loops, r_major: f64, r: f64, current: f64, samples: usize) -> (f64, [f64; 3]) {
    let pts: Vec<([f64; 3], [f64; 3])> =
        (0..samples).map(|j| loop_point(r_major, r, 2.0 * PI * (j as f64 + 0.5) / samples as f64)).collect();
    let b = loops.b(&pts.iter().map(|p| p.0).collect::<Vec<_>>());
    let ds = 2.0 * PI * r / samples as f64;
    let (mut radial, mut f) = (0.0, [0.0; 3]);
    for ((p, t), b) in pts.iter().zip(&b) {
        let df = cross(*t, *b);
        let rho_hat = [(p[0] - r_major) / r, 0.0, p[2] / r];
        radial += current * dot(df, rho_hat) * ds;
        for i in 0..3 {
            f[i] += current * df[i] * ds;
        }
    }
    (radial, f)
}

pub fn run() -> bool {
    let mut t = Tally { pass: 0, fail: 0 };
    let sh = design::baseline();
    let loops = sh.loops();
    let i = sh.current();
    println!(
        "Baseline: tube centre {:.0} m, habitat {:.0} m, windings {:.0} / {:.0} m, {} loops each, {:.3e} A-turns ({:.3e} A a loop)",
        sh.r_major, design::R_HAB, sh.r1, sh.r2, sh.n, sh.amp_turns, i
    );
    println!("          field {:.2} T at the tube centre's radius, bending power {:.1} T m\n", sh.b_at(sh.r_major), sh.bending_power());

    println!("cfsem field of the windings:");
    // 1. One loop: the field at its centre is mu0 I / 2a along +normal.
    let mut one = Loops::default();
    one.push(10.0, 1000.0, [3.0, -2.0, 1.0], [0.0, 1.0, 0.0], 0.0);
    let b = one.b(&[[3.0, -2.0, 1.0]])[0];
    let want = MU_0 * 1000.0 / (2.0 * 10.0);
    t.err("one loop: centre field mu0 I / 2a, +normal", &format!("{:.4e} T", b[1]), (b[1] - want).abs() / want + b[0].abs() + b[2].abs(), 1e-9);

    // 2. Discrete windings against the ideal field mid-shell, at a coil and between coils.
    let mid = 0.5 * (sh.r1 + sh.r2);
    let half = PI / sh.n as f64;
    let mut worst: f64 = 0.0;
    let mut pts = Vec::new();
    for th in [0.0, 0.5 * PI, PI, 1.5 * PI] {
        for ph in [0.0, half] {
            let (rr, z) = (sh.r_major + mid * th.cos(), mid * th.sin());
            pts.push([rr * ph.cos(), rr * ph.sin(), z]);
        }
    }
    for (p, b) in pts.iter().zip(loops.b(&pts)) {
        let bi = sh.b_ideal(*p);
        let e = [b[0] - bi[0], b[1] - bi[1], b[2] - bi[2]];
        worst = worst.max(norm(e) / norm(bi));
    }
    t.err("mid-shell field vs mu0 N I / 2 pi R (8 points)", "", worst, 1e-3);

    // 3. Inside the habitat the windings cancel; next to the inner winding some leaks.
    let centre = loops.b(&[[sh.r_major, 0.0, 0.0]])[0];
    let bc = sh.b_at(sh.r_major);
    t.err("field at the habitat's centre / shell field", &format!("{:.2e} T", norm(centre)), norm(centre) / bc, 1e-6);
    // Leakage at the habitat's surface: a row of currents spaced s leaks ~ exp(-2 pi x / s) of
    // the sheet's field at distance x, so finer windings close it off fast.
    for n in [1000usize, 2000, 4000] {
        let s = Shield { n, ..sh.clone() };
        let l = s.loops();
        let h = PI / n as f64;
        let mut lpts = Vec::new();
        for th in [0.0, 0.5 * PI, PI] {
            for ph in [0.0, h] {
                let (rr, z) = (s.r_major + design::R_HAB * th.cos(), design::R_HAB * th.sin());
                lpts.push([rr * ph.cos(), rr * ph.sin(), z]);
            }
        }
        let leak = l.b(&lpts).into_iter().map(norm).fold(0.0, f64::max);
        let s_in = 2.0 * PI * (s.r_major - s.r1) / n as f64;
        let est = s.b_at(s.r_major - s.r1) * (-2.0 * PI * (s.r1 - design::R_HAB) / s_in).exp();
        println!(
            "  info leakage at the habitat surface, {n} loops ({:.1} m apart inside): up to {:.3e} mT (exp(-2 pi x / s) estimate {:.3e} mT)",
            s_in,
            leak * 1e3,
            est * 1e3
        );
    }

    // 4. Stored energy: (1/2) sum of I x flux, against the closed form.
    let flux_in = loop_flux(&loops, sh.r_major, sh.r1, 4096);
    let flux_out = loop_flux(&loops, sh.r_major, sh.r2, 4096);
    let e = 0.5 * sh.n as f64 * (-i * flux_in + i * flux_out);
    let ei = sh.energy_ideal();
    t.err("stored energy (1/2) sum I Phi vs closed form", &format!("{:.4e} J vs {:.4e}", e, ei), (e - ei).abs() / ei, 2e-2);

    // 5. Forces: magnetic pressure B^2 / 2 mu0 pushes the inner winding in and the outer out;
    //    the 1/R field pulls each coil toward the spin axis by (1/N) dE/dR.
    let ni = sh.amp_turns;
    let nf = sh.n as f64;
    for (r, cur, sign, name) in [(sh.r1, -i, -1.0, "inner"), (sh.r2, i, 1.0, "outer")] {
        let (radial, _) = loop_forces(&loops, sh.r_major, r, cur, 4096);
        let want = sign * MU_0 * ni * ni * r / (2.0 * nf * (sh.r_major.powi(2) - r * r).sqrt());
        t.err(&format!("{name} coil: radial force vs pressure B^2/2mu0"), &format!("{:.4e} N", radial), (radial - want).abs() / want.abs(), 2e-2);
    }
    // The pair's net pull is a small difference of large forces, so its discrete error falls
    // with the loop spacing: check it converges as 1/N toward (1/N) dE/dR (scaled per loop).
    let r = sh.r_major;
    let mut errs = Vec::new();
    for n in [1000usize, 2000, 4000] {
        let s = Shield { n, ..sh.clone() };
        let l = s.loops();
        let cur = s.current();
        let (_, fi) = loop_forces(&l, r, s.r1, -cur, 8192);
        let (_, fo) = loop_forces(&l, r, s.r2, cur, 8192);
        let want = MU_0 * ni * ni / (2.0 * n as f64) * r * (1.0 / (r * r - s.r1 * s.r1).sqrt() - 1.0 / (r * r - s.r2 * s.r2).sqrt());
        let e = (fi[0] + fo[0] - want).abs() / want.abs();
        println!("  info centering force, {n} loops: {:.4e} N vs {:.4e} (err {:.2e})", fi[0] + fo[0], want, e);
        errs.push(e);
    }
    let order = (errs[0] / errs[2]).log2() / 2.0;
    t.err("centering force converges ~1/N to (1/N) dE/dR", &format!("order {order:.2}"), (order - 1.0).abs(), 0.25);
    t.err("centering force at 4000 loops", "", errs[2], 3e-2);

    println!("\nParticle tracer:");
    // 6. Uniform field: a proton of 1 GV turns on r = B rho / B; the Boris rotation keeps speed.
    let bz = 2.0;
    let br = brho(1e9);
    let (mut p, mut d) = ([0.0; 3], [1.0, 0.0, 0.0]);
    let r_l = br / bz;
    let steps = 20_000;
    let ds = 2.0 * PI * r_l / steps as f64;
    let (mut ymax, mut speed_err): (f64, f64) = (0.0, 0.0);
    for _ in 0..steps {
        d = boris_turn(d, [0.0, 0.0, bz], ds, br);
        p = [p[0] + ds * d[0], p[1] + ds * d[1], p[2] + ds * d[2]];
        ymax = ymax.max(p[1].abs());
        speed_err = speed_err.max((norm(d) - 1.0).abs());
    }
    t.err("gyro-diameter 2 B rho / B", &format!("{:.4} m vs {:.4}", ymax, 2.0 * r_l), (ymax - 2.0 * r_l).abs() / (2.0 * r_l), 1e-3);
    t.err("closes after one turn (2 pi B rho / B of path)", &format!("{:.2e} m", norm(p)), norm(p) / r_l, 1e-3);
    t.err("speed kept (|direction| = 1)", "", speed_err, 1e-12);
    println!("  info 1 GV proton: B rho {:.4} T m (c = {C} m/s)", br);

    // 7. A shell far from the spin axis and wide round the tube is a flat slab: at normal incidence a particle gets
    //    through only if its gyroradius exceeds the thickness, R > c B d.
    let slab = Shield::with_field(1.0e9, 1.0e6, 1.0e6 + 20.0, 1, 5.0);
    let tr = Tracer { r_hab: 1.0e6 - 5.0, ..Tracer::default() };
    let cut = C * slab.b_at(slab.r_major + 1.0e6 + 10.0) * 20.0; // V
    let start = [slab.r_major + tr.r_hab + 1e-3, 0.0, 0.0];
    let mut ok = 0;
    for (f, want) in [(0.97, Fate::Habitat), (1.03, Fate::Escaped)] {
        let fate = tr.follow(&slab, start, [1.0, 0.0, 0.0], -brho(f * cut));
        if fate == want {
            ok += 1;
        }
        println!("  info slab: {:.2} x c B d = {:.2} GV -> {:?} (want {:?})", f, f * cut / 1e9, fate, want);
    }
    t.err("slab cutoff at c B d (normal incidence)", &format!("{:.1} GV", cut / 1e9), (2 - ok) as f64, 0.0);

    println!("\nCosmic-ray dose model:");
    // 8. Bethe in water: a 100 MeV proton loses 7.289 MeV cm^2/g (NIST PSTAR); the minimum
    //    ionisation is 1.981 MeV cm^2/g (PDG water, with its density correction).
    let p = &dose::SPECIES[0];
    let s100 = dose::stopping(p, 0.1);
    t.err("proton stopping at 100 MeV vs PSTAR 7.289", &format!("{s100:.3} MeV cm^2/g"), (s100 - 7.289).abs() / 7.289, 1e-2);
    let smin = (0..400).map(|k| dose::stopping(p, 0.1 * 1.03f64.powi(k))).fold(f64::INFINITY, f64::min);
    t.err("proton minimum ionisation vs PDG 1.981", &format!("{smin:.3} MeV cm^2/g"), (smin - 1.981).abs() / 1.981, 3e-2);
    // 9. Free space at solar minimum: Ambroglini et al. 2016 (CREME 2009, Geant3) give a skin dose
    //    equivalent of 117.7 cSv/yr.
    let one = |_: f64| 1.0;
    let (hmin, dmin) = dose::total(0.42, &one);
    let (hmax, _) = dose::total(1.10, &one);
    t.err("free-space skin dose equivalent vs 117.7 cSv/yr", &format!("{:.1} cSv/yr (mean Q {:.2})", hmin * 100.0, hmin / dmin), (hmin - 1.177).abs() / 1.177, 0.35);
    // Skin and organ doses modulate differently (low-energy heavy ions dominate the skin and
    // stop in the body), so the solar cycle here is information only.
    println!("  info solar max (1.1 GV) / min (0.42 GV) skin dose equivalent {:.3} (organ dose, Ambroglini: 39.5 / 62.4 = 0.633)", hmax / hmin);
    for s in dose::SPECIES.iter() {
        let (h, _) = dose::species_dose(s, 0.42, &one);
        println!("  info {:<5} {:.1} cSv/yr", s.name, h * 100.0);
    }

    println!("\n{} passed, {} failed", t.pass, t.fail);
    t.fail == 0
}
