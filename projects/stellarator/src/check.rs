//! `check`: the coil model, cfsem's straight-filament field for these coils, and the field-line
//! tracer against known results before any design is built on them.

use std::f64::consts::PI;
use std::time::Instant;

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_density_circular_filament_cartesian_scalar;

use crate::coils::{Filaments, RotatingEllipse};
use crate::design;
use crate::trace::Tracer;

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
        println!("  {tag} {what:<46} {detail} err {err:.2e} (tol {tol:.0e})");
    }
}

/// Ellipse perimeter by the Gauss-Kummer series in h = ((a - b) / (a + b))^2.
fn ellipse_perimeter(a: f64, b: f64) -> f64 {
    let h = ((a - b) / (a + b)).powi(2);
    let (mut sum, mut binom, mut hn) = (1.0, 1.0, 1.0);
    for n in 1..60 {
        // binom(1/2, n) from binom(1/2, n - 1).
        binom *= (0.5 - (n as f64 - 1.0)) / n as f64;
        hn *= h;
        sum += binom * binom * hn;
    }
    PI * (a + b) * sum
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

pub fn run() -> bool {
    let mut t = Tally { pass: 0, fail: 0 };
    let base = design::baseline();

    println!("Fourier coils");
    {
        let c = base.coil(1);
        let n = 64;
        let (a, b) = (base.a_long, base.a_short);
        let l = c.length(n);
        let p = ellipse_perimeter(a, b);
        t.err(
            "order-1 coil length = ellipse perimeter",
            &format!("{l:.12} m vs {p:.12} m"),
            (l / p - 1.0).abs(),
            1e-12,
        );
        // Every sampled point satisfies the ellipse equation in its own (u, v) frame.
        let centre = c.point(0.0);
        let opposite = c.point(PI);
        let mid: [f64; 3] = std::array::from_fn(|i| 0.5 * (centre[i] + opposite[i]));
        let u: [f64; 3] = std::array::from_fn(|i| (centre[i] - mid[i]) / a);
        let q = c.point(0.5 * PI);
        let v: [f64; 3] = std::array::from_fn(|i| -(q[i] - mid[i]) / b);
        let mut e = 0.0_f64;
        for p in c.sample(97) {
            let d: [f64; 3] = std::array::from_fn(|i| p[i] - mid[i]);
            let pu = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
            let pv = d[0] * v[0] + d[1] * v[1] + d[2] * v[2];
            let off = norm(std::array::from_fn(|i| d[i] - pu * u[i] - pv * v[i]));
            e = e.max(((pu / a).powi(2) + (pv / b).powi(2) - 1.0).abs()).max(off / a);
        }
        t.err("points on the ellipse, in its plane", "97 points", e, 1e-14);
        // Tangent against a central difference of the position.
        let (tt, h) = (0.7, 1e-5);
        let (p1, p0, d) = (c.point(tt + h), c.point(tt - h), c.tangent(tt));
        let fd: [f64; 3] = std::array::from_fn(|i| (p1[i] - p0[i]) / (2.0 * h));
        let e = norm(std::array::from_fn(|i| fd[i] - d[i])) / norm(d);
        t.err("tangent = d point / dt", "central difference", e, 1e-9);
    }

    println!("cfsem filaments vs circular filament (one circular coil)");
    {
        // A circular coil (a_long = a_short, no shift): the polygon of linear filaments
        // against cfsem's oriented circular filament, at points inside and outside the coil.
        let circ = RotatingEllipse {
            a_short: base.a_long,
            shift: 0.0,
            shift_z: 0.0,
            ..base
        };
        let k = 0;
        let coil = circ.coil(k);
        let phi = circ.phi(k);
        let (s, c) = phi.sin_cos();
        let pts = [
            [circ.rc * c, circ.rc * s, 0.0],
            [(circ.rc + 0.05) * c - 0.02 * s, (circ.rc + 0.05) * s + 0.02 * c, 0.03],
            [0.4 * c + 0.1 * s, 0.4 * s - 0.1 * c, -0.1],
        ];
        let mut errs = Vec::new();
        for segs in [64, 256, 1024] {
            let f = Filaments::new(std::slice::from_ref(&coil), circ.current, segs, 0.0);
            let mut e = 0.0_f64;
            for p in pts {
                let got = f.b(p);
                let w = flux_density_circular_filament_cartesian_scalar(
                    (circ.a_long, circ.current),
                    (circ.rc * c, circ.rc * s, 0.0),
                    (-s, c, 0.0),
                    0.0,
                    (p[0], p[1], p[2]),
                );
                let want = [w.0, w.1, w.2];
                e = e.max(norm(std::array::from_fn(|i| got[i] - want[i])) / norm(want));
            }
            errs.push(e);
            println!("       {segs:>5} segments: max rel error {e:.2e}");
        }
        let order = (errs[0] / errs[2]).log2() / 4.0;
        t.err("polygon converges at order 2", &format!("order {order:.3}"), (order - 2.0).abs(), 0.05);
        t.err("1024 segments", "", errs[2], 1e-5);
    }

    let f = design::filaments(&base);
    println!(
        "baseline coils: {} coils x {} segments = {} filaments, {:.1} kA-turns each",
        base.n_coils(),
        design::SEGS,
        f.len(),
        base.current / 1e3
    );

    println!("Ampere's law");
    {
        // A toroidal circle links every coil once: integral B.dl = mu0 N I.
        let n = 720;
        for (r, z) in [(base.rc, 0.0), (base.rc - 0.01, 0.005)] {
            let mut sum = 0.0;
            for k in 0..n {
                let phi = 2.0 * PI * k as f64 / n as f64;
                sum += f.b_cyl(r, phi, z)[1] * r * 2.0 * PI / n as f64;
            }
            let want = MU_0 * base.n_coils() as f64 * base.current;
            t.err(
                "toroidal loop: mu0 N I",
                &format!("R {r:.3} Z {z:+.3}: {sum:.6e} vs {want:.6e} T m"),
                (sum / want - 1.0).abs(),
                1e-10,
            );
        }
        // A poloidal circle in the phi = 0 plane links no coil.
        let mut sum = 0.0;
        let rho = 0.04;
        for k in 0..n {
            let th = 2.0 * PI * k as f64 / n as f64;
            let b = f.b_cyl(base.rc + rho * th.cos(), 0.0, rho * th.sin());
            sum += (-b[0] * th.sin() + b[2] * th.cos()) * rho * 2.0 * PI / n as f64;
        }
        let scale = design::B0 * 2.0 * PI * rho;
        t.err("poloidal loop links nothing", &format!("{sum:.2e} T m"), (sum / scale).abs(), 1e-12);
    }

    println!("div B, curl B, symmetry");
    {
        let h = 1e-5;
        let pts = [
            [base.rc, 0.0, 0.0],
            [base.rc + 0.03, 0.2, 0.02],
            [-0.2, 0.15, -0.03],
        ];
        // Both with the conductor radius (cfsem blends inside it) and as bare filaments.
        let f0 = Filaments::new(&base.coils(), base.current, design::SEGS, 0.0);
        for (label, ff) in [("bare filaments", &f0), ("4 mm wire radius", &f)] {
        let (mut div_e, mut curl_e, mut worst) = (0.0_f64, 0.0_f64, [0.0; 3]);
        for p in pts {
            let d = |i: usize, s: f64| {
                let mut q = p;
                q[i] += s;
                ff.b(q)
            };
            let g: [[f64; 3]; 3] = std::array::from_fn(|i| {
                let (bp, bm) = (d(i, h), d(i, -h));
                std::array::from_fn(|j| (bp[j] - bm[j]) / (2.0 * h))
            });
            let scale = norm(ff.b(p)) / 0.05;
            div_e = div_e.max((g[0][0] + g[1][1] + g[2][2]).abs() / scale);
            let curl = [g[1][2] - g[2][1], g[2][0] - g[0][2], g[0][1] - g[1][0]];
            if norm(curl) / scale > curl_e {
                worst = p;
            }
            curl_e = curl_e.max(norm(curl) / scale);
        }
        println!("       {label}: worst curl at ({:.3}, {:.3}, {:.3})", worst[0], worst[1], worst[2]);
        t.err(&format!("div B = 0, {label}"), "relative to |B| / 5 cm", div_e, 1e-6);
        t.err(&format!("curl B = 0, {label}"), "relative to |B| / 5 cm", curl_e, 1e-6);
        }
        let per = 2.0 * PI / base.n_fp as f64;
        let (mut e_per, mut e_sym) = (0.0_f64, 0.0_f64);
        for (r, phi, z) in [(base.rc, 0.3, 0.01), (base.rc + 0.04, 1.1, -0.03), (base.rc - 0.02, 2.0, 0.05)] {
            let b = f.b_cyl(r, phi, z);
            let bp = f.b_cyl(r, phi + per, z);
            let bs = f.b_cyl(r, -phi, -z);
            let n = norm(b);
            e_per = e_per.max(norm(std::array::from_fn(|i| bp[i] - b[i])) / n);
            e_sym = e_sym.max(norm([bs[0] + b[0], bs[1] - b[1], bs[2] - b[2]]) / n);
        }
        t.err("field-period symmetry", "B(phi + 2pi/n_fp) = B(phi)", e_per, 1e-12);
        t.err("stellarator symmetry", "B_R odd, B_phi B_Z even", e_sym, 1e-12);
    }

    println!("field-line tracer");
    {
        let tr = design::tracer(&base);
        // RK4 order from the period map at 60 / 120 / 240 steps against 1920.
        let y0 = (base.rc + 0.03, 0.01);
        let reference = Tracer { steps: 1920, ..tr }.period_map(&f, y0);
        let mut errs = Vec::new();
        for steps in [60, 120, 240] {
            let y = Tracer { steps, ..tr }.period_map(&f, y0);
            errs.push((y.0 - reference.0).hypot(y.1 - reference.1));
        }
        let order = (errs[0] / errs[2]).log2() / 2.0;
        println!("       60 / 120 / 240 steps: {:.2e} / {:.2e} / {:.2e} m", errs[0], errs[1], errs[2]);
        // Over a whole period the leading error terms partly cancel, so the observed order runs above 4.
        t.err("RK4 converges at order >= 4", &format!("order {order:.3}"), (4.0 - order).max(0.0), 0.2);
        t.err(
            &format!("{} steps a period", tr.steps),
            "",
            {
                let y = tr.period_map(&f, y0);
                (y.0 - reference.0).hypot(y.1 - reference.1)
            },
            1e-7,
        );

        // Circular planar coils: an up-down symmetric axisymmetric-like field with ripple only.
        // The axis sits on Z = 0 and there is no rotational transform.
        let circ = RotatingEllipse {
            a_short: base.a_long,
            shift: 0.0,
            shift_z: 0.0,
            ..base
        };
        let fc = design::filaments(&circ);
        let tc = design::tracer(&circ);
        let (axis, res) = tc.find_axis(&fc, (circ.rc, 0.0));
        // Each polygon starts at the ellipse angle, so the polygons are up-down symmetric only to the
        // 128-segment error (~5e-4 of the field): the axis sits on Z = 0 to that order.
        t.err("circular coils: axis near Z = 0", &format!("R {:.6} m, residual {res:.1e}", axis.0), axis.1.abs() / circ.a_long, 5e-4);
        let path = tc.axis_path(&fc, axis);
        let line = tc.trace(&fc, (axis.0 + 0.03, 0.0), &path, 30);
        let (iota, _) = line.iota(circ.n_fp);
        t.err("circular coils: no rotational transform", &format!("iota {iota:.2e}"), iota.abs(), 1e-9);

        // The baseline axis: the elliptic fixed point nearest the coil centre. Stellarator
        // symmetry puts it on Z = 0 at phi = 0; RK4 is not time-symmetric, so the discrete map
        // keeps that symmetry only to the step error.
        let t0 = Instant::now();
        let ax = tr.axis(&f, (base.rc, 0.0)).expect("an elliptic axis");
        t.err(
            "baseline: axis on Z = 0 at phi = 0",
            &format!("R {:.6} m, residual {:.1e}, {:.2} s", ax.rz.0, ax.residual, t0.elapsed().as_secs_f64()),
            ax.rz.1.abs(),
            1e-6,
        );
        // A field-line map preserves the flux B_phi dR dZ; at a fixed point that makes det J = 1.
        // RK4 does not preserve it exactly: the error should fall with the step.
        let j2 = Tracer { steps: 2 * tr.steps, ..tr }.jacobian(&f, ax.rz);
        let det2 = j2[0][0] * j2[1][1] - j2[0][1] * j2[1][0];
        println!("       det J - 1: {:.2e} at {} steps, {:.2e} at {}", ax.det - 1.0, tr.steps, det2 - 1.0, 2 * tr.steps);
        t.err("period map preserves flux: det J = 1", &format!("{} steps", 2 * tr.steps), (det2 - 1.0).abs(), 2e-7);
        // Rotational transform on the axis two ways: the Jacobian's eigen-rotation, and a line
        // 2 mm off the axis traced for 600 periods (Birkhoff average).
        let path = tr.axis_path(&f, ax.rz);
        let line = tr.trace(&f, (ax.rz.0 + 0.002, 0.0), &path, 600);
        let (iota, mean) = line.iota(base.n_fp);
        t.err(
            "iota near the axis: traced vs Jacobian",
            &format!("{iota:.6} (mean {mean:.6}) vs {:.6}", ax.iota(base.n_fp)),
            (iota.abs() / ax.iota(base.n_fp) - 1.0).abs(),
            2e-3,
        );
    }

    println!("{} passed, {} failed", t.pass, t.fail);
    t.fail == 0
}
