//! `check`: cfsem's straight-filament inductance and field against textbook results, before any
//! design is built on them.

use std::f64::consts::PI;

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_circular_filament;

use crate::loops::{Loop, P3, field, mutual, self_inductance, self_inductance_raw};

/// (H) partial self-inductance of a straight round wire, length l, radius a: external (flux
/// outside the wire, exact for a thin straight wire) and the uniform-current internal term.
fn partial_self(l: f64, a: f64) -> (f64, f64) {
    let s = (l * l + a * a).sqrt();
    let ext = MU_0 / (2.0 * PI) * (l * ((l + s) / a).ln() - s + a);
    (ext, MU_0 / (8.0 * PI) * l)
}

/// (H) partial mutual inductance of two parallel, side-by-side filaments of length l, distance d.
fn partial_mutual(l: f64, d: f64) -> f64 {
    let s = (l * l + d * d).sqrt();
    MU_0 / (2.0 * PI) * (l * ((l + s) / d).ln() - s + d)
}

/// (H) a rectangle of round wire (Grover / Ruehli partial inductances; perpendicular sides do
/// not couple): (external only, with the DC internal term).
fn rect_closed(lx: f64, ly: f64, a: f64) -> (f64, f64) {
    let (ex, ix) = partial_self(lx, a);
    let (ey, iy) = partial_self(ly, a);
    let ext = 2.0 * (ex + ey) - 2.0 * partial_mutual(lx, ly) - 2.0 * partial_mutual(ly, lx);
    (ext, ext + 2.0 * (ix + iy))
}

/// (H) Neumann double integral between two closed polylines, Gauss-Legendre on both, written
/// independently of cfsem (filaments; no wire radius, so for separate loops only).
fn neumann(a: &Loop, b: &Loop, q: usize) -> f64 {
    // Gauss-Legendre nodes and weights on [0, 1] by Newton on P_q.
    let mut gl = Vec::new();
    for k in 0..q {
        let mut x = (PI * (k as f64 + 0.75) / (q as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for n in 2..=q {
                let p2 = ((2 * n - 1) as f64 * x * p1 - (n - 1) as f64 * p0) / n as f64;
                p0 = p1;
                p1 = p2;
            }
            dp = q as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        gl.push((0.5 * (1.0 + x), 1.0 / ((1.0 - x * x) * dp * dp)));
    }
    let (sa, da) = a.segments();
    let (sb, db) = b.segments();
    let mut m = 0.0;
    for (pa, va) in sa.iter().zip(&da) {
        for (pb, vb) in sb.iter().zip(&db) {
            let dot = va.0 * vb.0 + va.1 * vb.1 + va.2 * vb.2;
            if dot == 0.0 {
                continue;
            }
            let mut s = 0.0;
            for &(ta, wa) in &gl {
                for &(tb, wb) in &gl {
                    let r = (
                        pa.0 + ta * va.0 - pb.0 - tb * vb.0,
                        pa.1 + ta * va.1 - pb.1 - tb * vb.1,
                        pa.2 + ta * va.2 - pb.2 - tb * vb.2,
                    );
                    s += wa * wb / (r.0 * r.0 + r.1 * r.1 + r.2 * r.2).sqrt();
                }
            }
            m += dot * s;
        }
    }
    MU_0 / (4.0 * PI) * a.turns * b.turns * m
}

pub fn run() -> bool {
    let mut ok = true;
    let mut check = |name: &str, got: f64, want: f64, tol: f64| {
        let e = got / want - 1.0;
        let pass = e.abs() <= tol;
        ok &= pass;
        println!(
            "   {} {name}: {got:.6e} vs {want:.6e} ({e:+.2e}, tol {tol:.0e})",
            if pass { "ok  " } else { "FAIL" }
        );
    };

    println!(
        "1. Self-inductance of a 100 x 50 mm rectangle of 1 mm wire (cfsem A.dl on the wire axis)"
    );
    let (lx, ly, a) = (0.1, 0.05, 0.0005);
    let (ext, dc) = rect_closed(lx, ly, a);
    println!(
        "   closed form (Grover partial inductances): external {:.4} nH, with the DC internal term {:.4} nH",
        ext * 1e9,
        dc * 1e9
    );
    for n in [1, 4, 16, 64, 256] {
        let rect = Loop::rect(0.0, 0.0, 0.0, lx, ly, n, 1.0);
        let l = self_inductance_raw(&rect, a);
        println!(
            "   {n:>3} segments a side: cfsem {:.4} nH; vs external {:+.2e}, vs DC {:+.2e}, internal implied {:.3} x mu0/8pi per metre",
            l * 1e9,
            l / ext - 1.0,
            l / dc - 1.0,
            (l - ext) / (MU_0 / (8.0 * PI) * rect.length())
        );
        if n == 64 {
            check(
                "corrected (less mu0/8pi per metre), 64 a side",
                self_inductance(&rect, a),
                dc,
                2e-3,
            );
        }
    }

    println!(
        "2. Self-inductance of a 50 mm radius ring of 1 mm wire, polygon vs Wien mu0 R (ln 8R/a - 7/4) (cfsem's own Python formula)"
    );
    let (r, a) = (0.05, 0.0005);
    let ra2 = (a / r) * (a / r);
    let wien = MU_0 * r * ((1.0 + 0.125 * ra2) * (8.0 * r / a).ln() - 0.0083 * ra2 - 1.75);
    for m in [16, 64, 256, 1024] {
        let ring = Loop::circle((0.0, 0.0, 0.0), r, m, 1.0);
        let l = self_inductance_raw(&ring, a);
        println!(
            "   {m:>4}-gon: cfsem {:.4} nH vs Wien {:.4} nH ({:+.2e}); excess {:.4} nH, mu0/8pi x length {:.4} nH",
            l * 1e9,
            wien * 1e9,
            l / wien - 1.0,
            (l - wien) * 1e9,
            MU_0 / (8.0 * PI) * ring.length() * 1e9
        );
        if m == 1024 {
            check("corrected, 1024-gon", self_inductance(&ring, a), wien, 1e-4);
        }
    }

    println!(
        "3. Mutual inductance of coaxial circles (R 50 and 40 mm, 20 mm apart): polygons vs cfsem's circular filament"
    );
    let mut psi = [0.0];
    flux_circular_filament(
        (&[0.04], &[0.02], &[1.0]),
        &[0.0],
        (&[0.05], &[0.0]),
        &mut psi,
    )
    .expect("lengths");
    for m in [16, 64, 256, 1024] {
        let got = mutual(
            &Loop::circle((0.0, 0.0, 0.0), 0.05, m, 1.0),
            &Loop::circle((0.0, 0.0, 0.02), 0.04, m, 1.0),
        );
        if m == 1024 {
            check("1024-gons", got, psi[0], 1e-5);
        } else {
            println!("   {m:>4}-gons: {got:.6e} ({:+.2e})", got / psi[0] - 1.0);
        }
    }

    println!(
        "4. Mutual inductance of two squares (60 mm and 40 mm, 15 mm apart, 10 mm offset): cfsem vs own Neumann quadrature"
    );
    let a0 = Loop::rect(0.0, 0.0, 0.0, 0.06, 0.06, 8, 3.0);
    let b0 = Loop::rect(0.01, 0.004, 0.015, 0.04, 0.04, 8, 2.0);
    let want = neumann(&a0, &b0, 24);
    check("squares", mutual(&a0, &b0), want, 1e-6);
    check("symmetric", mutual(&b0, &a0), want, 1e-6);

    println!(
        "5. Field of a 20 m track loop (conductors 50 mm apart) at mid-length: finite wires (Biot-Savart closed form), and the 2D two-wire limit"
    );
    let track = Loop::rect(0.0, 0.0, 0.0, 20.0, 0.05, 1, 1.0);
    for p in [(0.0, 0.0, 0.02), (0.0, 0.04, 0.01), (0.0, -0.1, 0.03)] {
        let b = field(&[(&track, 1.0)], &[p])[0];
        let (mut want, mut by2, mut bz2) = ((0.0, 0.0, 0.0), 0.0, 0.0);
        let (s, d) = track.segments();
        for (a0, dl) in s.iter().zip(&d) {
            // Straight segment: B = mu0 I / (4 pi rho) (cos th1 - cos th2) along dl x rho_hat.
            let len = (dl.0 * dl.0 + dl.1 * dl.1 + dl.2 * dl.2).sqrt();
            let u = (dl.0 / len, dl.1 / len, dl.2 / len);
            let w = (p.0 - a0.0, p.1 - a0.1, p.2 - a0.2);
            let t = w.0 * u.0 + w.1 * u.1 + w.2 * u.2;
            let rho = (w.0 - t * u.0, w.1 - t * u.1, w.2 - t * u.2);
            let rr = (rho.0 * rho.0 + rho.1 * rho.1 + rho.2 * rho.2).sqrt();
            let k =
                MU_0 / (4.0 * PI * rr) * (t / t.hypot(rr) - (t - len) / (t - len).hypot(rr)) / rr;
            want.0 += k * (u.1 * rho.2 - u.2 * rho.1);
            want.1 += k * (u.2 * rho.0 - u.0 * rho.2);
            want.2 += k * (u.0 * rho.1 - u.1 * rho.0);
        }
        for (y0, i) in [(-0.025, 1.0), (0.025, -1.0)] {
            let (dy, dz): (f64, f64) = (p.1 - y0, p.2);
            let k = MU_0 * i / (2.0 * PI * (dy * dy + dz * dz));
            by2 += -k * dz;
            bz2 += k * dy;
        }
        let err = (b.0 - want.0).hypot(b.1 - want.1).hypot(b.2 - want.2) / want.1.hypot(want.2);
        check(
            &format!("|B| at y {:+.0} mm, z {:.0} mm", p.1 * 1e3, p.2 * 1e3),
            b.1.hypot(b.2),
            want.1.hypot(want.2),
            1e-9,
        );
        println!(
            "        vector error {err:.1e}; the 2D infinite-wire value differs by {:+.1e} (the 10 m ends)",
            by2.hypot(bz2) / want.1.hypot(want.2) - 1.0
        );
    }

    println!(
        "6. Images across z = 0 (a 30 mm square at z = 10 mm): boundary conditions on the plane"
    );
    let src = Loop::rect(0.003, -0.002, 0.01, 0.03, 0.03, 4, 1.0);
    let probe: Vec<P3> = (0..40)
        .map(|k| (-0.03 + 0.0015 * k as f64, 0.004 - 0.0004 * k as f64, 0.0))
        .collect();
    let free = field(&[(&src, 1.0)], &probe);
    let scale = free
        .iter()
        .map(|b| (b.0 * b.0 + b.1 * b.1 + b.2 * b.2).sqrt())
        .fold(0.0, f64::max);
    let cond = field(&[(&src, 1.0), (&src.image(), -1.0)], &probe);
    let iron = field(&[(&src, 1.0), (&src.image(), 1.0)], &probe);
    let worst_n = cond.iter().map(|b| b.2.abs()).fold(0.0, f64::max) / scale;
    let worst_t = iron.iter().map(|b| b.0.hypot(b.1)).fold(0.0, f64::max) / scale;
    let pass = worst_n < 1e-12 && worst_t < 1e-12;
    ok &= pass;
    println!(
        "   {} perfect conductor: normal B on the plane {worst_n:.1e} of the peak; ideal iron: tangential B {worst_t:.1e}",
        if pass { "ok  " } else { "FAIL" }
    );
    ok
}
