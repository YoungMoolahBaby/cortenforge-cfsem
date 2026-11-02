//! `check`: cfsem's Grad-Shafranov operator, and the flux and field conventions, against known
//! results before any equilibrium is built on them.

use std::time::Instant;

use cfsem::physics::circular_filament::{flux_circular_filament, flux_density_circular_filament};

use crate::gs::{Grid, Solver, apply, field_from_flux, operator, rhs};

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
        println!("  {tag} {what:<44} {detail} err {err:.2e} (tol {tol:.0e})");
    }
}

/// Max |numeric - exact| over interior nodes at least `margin` from the edge, relative to
/// max |exact| over the same nodes.
fn rel_max(grid: &Grid, got: &[f64], want: &[f64], margin: usize) -> f64 {
    let (mut e, mut m) = (0.0_f64, 0.0_f64);
    for i in margin..grid.nr() - margin {
        for j in margin..grid.nz() - margin {
            let k = grid.idx(i, j);
            e = e.max((got[k] - want[k]).abs());
            m = m.max(want[k].abs());
        }
    }
    e / m
}

/// Test functions with exact Delta* = f_RR - f_R / R + f_ZZ.
struct TestFn {
    name: &'static str,
    f: fn(f64, f64) -> f64,
    lap: fn(f64, f64) -> f64,
}

const TESTS: [TestFn; 3] = [
    TestFn {
        name: "R^4 (Delta* = 8 R^2)",
        f: |r, _| r.powi(4),
        lap: |r, _| 8.0 * r * r,
    },
    TestFn {
        name: "R^2 Z^2 (Delta* = 2 R^2)",
        f: |r, z| r * r * z * z,
        lap: |r, _| 2.0 * r * r,
    },
    TestFn {
        name: "R^2 cos R exp(-Z^2)",
        f: |r, z| r * r * r.cos() * (-z * z).exp(),
        lap: |r, z| {
            let g = (-z * z).exp();
            g * (-3.0 * r * r.sin() - r * r * r.cos() + r * r * r.cos() * (4.0 * z * z - 2.0))
        },
    },
];

pub fn run() -> bool {
    let mut t = Tally { pass: 0, fail: 0 };

    println!(
        "1. Delta* applied to known functions, R 0.5-1.5, Z -0.5-0.5 (2 nodes in from the edge)"
    );
    for tf in &TESTS {
        for order4 in [false, true] {
            let mut errs = Vec::new();
            for n in [33, 65] {
                let g = Grid::new((0.5, 1.5), n, (-0.5, 0.5), n);
                let got = apply(&operator(&g, order4), &g.sample(tf.f));
                errs.push(rel_max(&g, &got, &g.sample(tf.lap), 2));
            }
            let observed = (errs[0] / errs[1]).log2();
            let label = format!("{} order {}", tf.name, if order4 { 4 } else { 2 });
            let detail = format!(
                "n 33 -> 65: {:.2e} -> {:.2e}, order {observed:.2}",
                errs[0], errs[1]
            );
            // These polynomials are exact for both stencils (the 2nd-order one is in flux form with
            // half-point radii), up to rounding that grows with n.
            if !tf.name.contains("cos") {
                t.err(&label, &detail, errs[1], 1e-9);
            } else {
                let want = if order4 { 4.0 } else { 2.0 };
                t.err(&label, &detail, (observed - want).abs() / want, 0.1);
            }
        }
    }

    println!(
        "2. Fixed-boundary Solov'ev equilibrium, Psi = R^2 Z^2 + (a/4)(R^2 - R0^2)^2, a = 1.5"
    );
    println!("   (Delta* Psi = 2 (1 + a) R^2; nested surfaces around the axis at (R0, 0), R0 = 1)");
    let (a, r0) = (1.5, 1.0);
    let sol = |r: f64, z: f64| r * r * z * z + 0.25 * a * (r * r - r0 * r0).powi(2);
    for order4 in [false, true] {
        let mut errs = Vec::new();
        for n in [33, 65] {
            let g = Grid::new((0.6, 1.4), n, (-0.5, 0.5), n);
            let b: Vec<f64> = g
                .nodes()
                .map(|(r, z, k)| {
                    let (i, j) = (k / g.nz(), k % g.nz());
                    if g.is_boundary(i, j) {
                        sol(r, z)
                    } else {
                        2.0 * (1.0 + a) * r * r
                    }
                })
                .collect();
            let psi = Solver::new(&g, order4).solve(&b);
            errs.push(rel_max(&g, &psi, &g.sample(sol), 1));
            if n == 65 {
                // The magnetic axis: the interior minimum of Psi.
                let k = (0..g.len())
                    .min_by(|&x, &y| psi[x].total_cmp(&psi[y]))
                    .expect("nodes");
                let (ra, za) = (g.rs[k / g.nz()], g.zs[k % g.nz()]);
                t.err(
                    &format!("  axis, order {}", if order4 { 4 } else { 2 }),
                    &format!("at R {ra:.4}, Z {za:+.4}"),
                    (ra - r0).hypot(za),
                    1e-9,
                );
            }
        }
        let observed = (errs[0] / errs[1]).log2();
        let label = format!("solve, order {}", if order4 { 4 } else { 2 });
        let detail = format!(
            "n 33 -> 65: {:.2e} -> {:.2e}, order {observed:.2}",
            errs[0], errs[1]
        );
        // A degree-4 polynomial: both stencils reproduce it to rounding.
        t.err(&label, &detail, errs[1], 1e-9);
    }

    println!(
        "3. Vacuum flux of a coil: grid solve vs cfsem's filament flux (same discrete source)"
    );
    // A coil block of current density J on the nodes inside it; the same current as filaments at
    // those nodes is the exact source for the grid, so only the stencil differs.
    let j_coil = 1e6; // A/m^2
    for (n, order4) in [(81, false), (161, false), (81, true), (161, true)] {
        let g = Grid::new((0.2, 2.0), n, (-1.0, 1.0), n);
        let inside = |r: f64, z: f64| (1.2..=1.4).contains(&r) && (0.3..=0.5).contains(&z);
        let jphi = g.sample(|r, z| if inside(r, z) { j_coil } else { 0.0 });
        let fil: Vec<(f64, f64)> = g
            .nodes()
            .filter(|n| inside(n.0, n.1))
            .map(|n| (n.0, n.1))
            .collect();
        let (rf, zf): (Vec<f64>, Vec<f64>) = fil.iter().copied().unzip();
        let i_each = vec![j_coil * g.dr() * g.dz(); fil.len()];
        let wire = vec![0.0; fil.len()];
        let exact = |r: &[f64], z: &[f64]| {
            let mut out = vec![0.0; r.len()];
            flux_circular_filament((&rf, &zf, &i_each), &wire, (r, z), &mut out).expect("lengths");
            out
        };
        let (rn, zn): (Vec<f64>, Vec<f64>) = g.nodes().map(|n| (n.0, n.1)).unzip();
        let want = exact(&rn, &zn);
        let t0 = Instant::now();
        let solver = Solver::new(&g, order4);
        let factor = t0.elapsed();
        let psi = solver.solve(&rhs(&g, &jphi, |r, z| exact(&[r], &[z])[0]));
        let nans = psi.iter().filter(|x| !x.is_finite()).count();
        if nans > 0 {
            println!("    {nans} of {} solved values are not finite", psi.len());
        }
        // Compare away from the coil (3 cells) and the edge.
        let (mut e, mut m) = (0.0_f64, 0.0_f64);
        for (r, z, k) in g.nodes() {
            let far = (r - 1.3).abs() > 0.1 + 3.0 * g.dr() || (z - 0.4).abs() > 0.1 + 3.0 * g.dz();
            let (i, jj) = (k / g.nz(), k % g.nz());
            if far && !g.is_boundary(i, jj) {
                e = e.max((psi[k] - want[k]).abs());
                // Normalize over the same nodes: the ideal filaments are singular at their own.
                m = m.max(want[k].abs());
            }
        }
        let label = format!("{n} x {n}, order {}", if order4 { 4 } else { 2 });
        let detail = format!(
            "{} filaments, LU {:.0} ms",
            fil.len(),
            factor.as_secs_f64() * 1e3
        );
        t.err(&label, &detail, e / m, if order4 { 1e-4 } else { 1e-2 });
    }

    println!("4. B from Psi: 4th-order differences of cfsem's flux vs cfsem's field, one loop");
    // The loop sits on a node, where the ideal flux is infinite: only nodes >= 0.1 m away count,
    // for errors and for normalization alike.
    let loop_ = (vec![1.0], vec![0.2], vec![1e5]);
    let far = |r: f64, z: f64| (r - 1.0).hypot(z - 0.2) > 0.1;
    let mut b_errs = Vec::new();
    let mut lap_err = 0.0;
    for (nr, nz) in [(91, 101), (181, 201)] {
        let g = Grid::new((0.2, 2.0), nr, (-1.0, 1.0), nz);
        let (rn, zn): (Vec<f64>, Vec<f64>) = g.nodes().map(|n| (n.0, n.1)).unzip();
        let mut psi = vec![0.0; g.len()];
        flux_circular_filament((&loop_.0, &loop_.1, &loop_.2), &[0.0], (&rn, &zn), &mut psi)
            .expect("lengths");
        let (br, bz) = field_from_flux(&g, &psi);
        let (mut wr, mut wz) = (vec![0.0; g.len()], vec![0.0; g.len()]);
        flux_density_circular_filament(
            (&loop_.0, &loop_.1, &loop_.2),
            &[0.0],
            (&rn, &zn),
            (&mut wr, &mut wz),
        )
        .expect("lengths");
        let (mut e, mut m) = (0.0_f64, 0.0_f64);
        for (r, z, k) in g.nodes() {
            if br[k].is_finite() && far(r, z) {
                e = e.max((br[k] - wr[k]).hypot(bz[k] - wz[k]));
                m = m.max(wr[k].hypot(wz[k]));
            }
        }
        b_errs.push(e / m);
        // Delta* of the exact vacuum flux is zero away from the source (truncation only).
        let lap = apply(&operator(&g, true), &psi);
        let (mut worst, mut scale) = (0.0_f64, 0.0_f64);
        for (r, z, k) in g.nodes() {
            let (i, j) = (k / g.nz(), k % g.nz());
            if i > 2 && j > 2 && i + 3 < g.nr() && j + 3 < g.nz() && far(r, z) {
                worst = worst.max(lap[k].abs());
                scale = scale.max(psi[k].abs() / (g.dr() * g.dr()));
            }
        }
        lap_err = worst / scale;
    }
    let observed = (b_errs[0] / b_errs[1]).log2();
    t.err(
        "|B| from Psi vs direct, >= 0.1 m from the loop",
        &format!(
            "dr 0.02 -> 0.01: {:.2e} -> {:.2e}, order {observed:.2}",
            b_errs[0], b_errs[1]
        ),
        (observed - 4.0).abs() / 4.0,
        0.1,
    );
    t.err(
        "Delta* of cfsem's vacuum flux, >= 0.1 m away",
        "relative to Psi/dr^2",
        lap_err,
        1e-6,
    );

    println!("\n{} passed, {} failed", t.pass, t.fail);
    t.fail == 0
}
