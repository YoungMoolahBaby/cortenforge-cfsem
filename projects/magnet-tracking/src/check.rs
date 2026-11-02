//! `check`: the forward models against textbook results before any tracking is built on them.

use std::f64::consts::PI;

use cfsem::MU_0;

use crate::magnet::{Magnet, V3, cel};

/// Complete elliptic integrals K(k), E(k) by the arithmetic-geometric mean (independent of cel).
fn k_e(k: f64) -> (f64, f64) {
    let (mut a, mut g) = (1.0, (1.0 - k * k).sqrt());
    let mut sum = 0.5 * k * k;
    let mut pow = 0.5;
    for _ in 0..40 {
        let (an, gn) = (0.5 * (a + g), (a * g).sqrt());
        let c = 0.5 * (a - g);
        pow *= 2.0;
        sum += pow * c * c;
        a = an;
        g = gn;
        if c.abs() < 1e-17 {
            break;
        }
    }
    let kk = PI / (2.0 * a);
    (kk, kk * (1.0 - sum))
}

pub fn run() -> bool {
    let mut ok = true;
    let mut check = |name: &str, err: f64, tol: f64| {
        let pass = err.abs() <= tol;
        ok &= pass;
        println!(
            "   {} {name}: {err:.1e} (tol {tol:.0e})",
            if pass { "ok  " } else { "FAIL" }
        );
    };

    println!(
        "1. cfsem's point dipole vs B = mu0 / 4 pi (3 (m.r) r / r^5 - m / r^3), 200 random points, a tilted moment"
    );
    let mag = Magnet {
        radius: 3e-3,
        length: 6e-3,
        br: 1.3,
    };
    let u = V3::new(0.3, -0.5, 0.8).normalize();
    let pos = V3::new(0.01, -0.02, 0.005);
    let mut seed = 12345u64;
    let mut rnd = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    let mut worst = 0.0_f64;
    for _ in 0..200 {
        let obs = V3::new(rnd(), rnd(), rnd()) * 0.2;
        let r = obs - pos;
        let m = u * mag.moment();
        let want =
            (r * (3.0 * m.dot(&r)) / r.norm().powi(5) - m / r.norm().powi(3)) * (MU_0 / (4.0 * PI));
        worst = worst.max((mag.dipole(&pos, &u, &obs) - want).norm() / want.norm());
    }
    check("worst relative difference", worst, 1e-13);

    println!(
        "2. Bulirsch's cel against K(k) and E(k) by the AGM: cel(kc, 1, 1, 1) = K, cel(kc, 1, 1, kc^2) = E"
    );
    let mut worst = 0.0_f64;
    for k in [0.1_f64, 0.5, 0.9, 0.99, 0.9999] {
        let kc = (1.0 - k * k).sqrt();
        let (kk, ee) = k_e(k);
        worst = worst
            .max((cel(kc, 1.0, 1.0, 1.0) / kk - 1.0).abs())
            .max((cel(kc, 1.0, 1.0, kc * kc) / ee - 1.0).abs());
    }
    check("worst relative difference (k 0.1 to 0.9999)", worst, 1e-13);

    println!(
        "3. Derby and Olbert on the axis vs the textbook B_z = Br/2 [(z+b)/sqrt(a^2+(z+b)^2) - (z-b)/sqrt(a^2+(z-b)^2)]"
    );
    let mut worst = 0.0_f64;
    for z in [0.0, 1e-3, 2.9e-3, 3.1e-3, 5e-3, 2e-2, 1e-1] {
        let (a, b) = (mag.radius, 0.5 * mag.length);
        let want = 0.5
            * mag.br
            * ((z + b) / (a * a + (z + b) * (z + b)).sqrt()
                - (z - b) / (a * a + (z - b) * (z - b)).sqrt());
        let (br, bz) = mag.derby(1e-12, z);
        worst = worst
            .max((bz / want - 1.0).abs())
            .max(br.abs() / want.abs());
    }
    check(
        "worst relative difference (z 0 to 100 mm, inside and out)",
        worst,
        1e-9,
    );

    println!(
        "4. cfsem current loops vs Derby and Olbert off the axis (tilted magnet, points 6-60 mm away)"
    );
    let pts: Vec<V3> = (0..60)
        .map(|k| {
            let t = k as f64 * 0.37;
            let r = 6e-3 + 54e-3 * (k as f64 / 59.0);
            pos + V3::new(
                t.cos() * (t * 0.7).sin(),
                t.sin() * (t * 0.7).sin(),
                (t * 0.7).cos(),
            ) * r
        })
        .collect();
    let exact: Vec<V3> = pts.iter().map(|p| mag.exact(&pos, &u, p)).collect();
    let mut last = 0.0;
    for n in [4, 16, 64, 256] {
        let b = mag.loops(&pos, &u, n, &pts);
        let worst = b
            .iter()
            .zip(&exact)
            .map(|(x, y)| (x - y).norm() / y.norm())
            .fold(0.0, f64::max);
        println!(
            "   {n:>3} loops: worst relative difference {worst:.2e}{}",
            if last > 0.0 {
                format!(" (x{:.1} smaller)", last / worst)
            } else {
                String::new()
            }
        );
        last = worst;
        if n == 256 {
            check("256 loops", worst, 1e-5);
        }
    }

    println!(
        "5. How far a point dipole is from the magnet, 45 deg off the axis (|B_dipole - B_exact| / |B_exact|)"
    );
    let along = V3::new(1.0, 0.0, 1.0).normalize();
    let z = V3::z();
    let o = V3::zeros();
    for d in [2.0, 4.0, 8.0, 16.0] {
        let p = along * d * 2.0 * mag.radius;
        let e =
            (mag.dipole(&o, &z, &p) - mag.exact(&o, &z, &p)).norm() / mag.exact(&o, &z, &p).norm();
        println!("   L/D 1.0, at {d:>2.0} diameters: {e:.2e}");
    }
    // The leading correction (the octupole) vanishes at L/D = sqrt(3) / 2 for a uniformly
    // magnetized cylinder; scan the aspect at fixed volume and distance.
    let vol = PI * 3e-3 * 3e-3 * 6e-3;
    let dist = 30e-3;
    let mut best = (0.0, f64::MAX);
    let mut row = Vec::new();
    for k in 0..=120 {
        let ld = 0.4 + 0.01 * k as f64;
        let r = (vol / (2.0 * PI * ld)).cbrt();
        let m = Magnet {
            radius: r,
            length: 2.0 * r * ld,
            br: 1.3,
        };
        let p = along * dist;
        let e = (m.dipole(&o, &z, &p) - m.exact(&o, &z, &p)).norm() / m.exact(&o, &z, &p).norm();
        if e < best.1 {
            best = (ld, e);
        }
        if k % 20 == 0 {
            row.push(format!("{ld:.1}: {e:.1e}"));
        }
    }
    println!("   same volume, 30 mm away: {}", row.join(", "));
    check(
        "aspect with the least dipole error vs sqrt(3)/2",
        best.0 - 3f64.sqrt() / 2.0,
        0.011,
    );
    println!("   (least error {:.1e} at L/D {:.2})", best.1, best.0);
    ok
}
