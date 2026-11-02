//! Reference formulas written here, independent of cfsem, to check it against.

use std::f64::consts::PI;

use cfsem::MU_0;

/// (Br, Bz) of an ideal loop (radius a, at z = 0, current i) by direct Biot-Savart quadrature
/// over the azimuth: trapezoid rule, which converges exponentially for this periodic integrand.
pub fn biot_savart_loop(a: f64, i: f64, (r, z): (f64, f64), n: usize) -> (f64, f64) {
    let (mut br, mut bz) = (0.0, 0.0);
    let dphi = 2.0 * PI / n as f64;
    for k in 0..n {
        let phi = k as f64 * dphi;
        let (c, s) = (phi.cos(), phi.sin());
        // Source a(cos, sin, 0), dl = a(-sin, cos, 0) dphi, R = obs - src with obs = (r, 0, z).
        let d = (r - a * c, -a * s, z);
        let d3 = (d.0 * d.0 + d.1 * d.1 + d.2 * d.2).powf(1.5);
        // dl x R = a (z cos, z sin, a - r cos) dphi; the y part integrates to zero.
        br += a * z * c / d3;
        bz += a * (a - r * c) / d3;
    }
    let k = MU_0 * i / (4.0 * PI) * dphi;
    (k * br, k * bz)
}

/// Complete elliptic integrals (K(m), E(m)) to machine precision by the arithmetic-geometric mean
/// (DLMF 19.8.1 and 19.8.6), independent of cfsem's curve fits.
pub fn ellip_agm(m: f64) -> (f64, f64) {
    ellip_agm_c(m, 1.0 - m)
}

/// [ellip_agm] given m and its complement 1 - m computed without cancellation, which keeps
/// K finite as m -> 1 (near-coincident rings).
pub fn ellip_agm_c(m: f64, mc: f64) -> (f64, f64) {
    let (mut a, mut b) = (1.0_f64, mc.sqrt());
    let mut sum = 0.5 * m; // 2^(n-1) c_n^2 summed from n = 0, c_0^2 = m
    let mut pow = 0.5;
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

/// (Br, Bz) of an ideal loop from the textbook elliptic-integral closed form (e.g. Simpson et al.
/// 2001, eqs. 3-4 in cylindrical form), evaluated with [ellip_agm].
pub fn closed_form_loop(a: f64, i: f64, (r, z): (f64, f64)) -> (f64, f64) {
    let p2 = (a + r) * (a + r) + z * z;
    let q2 = (a - r) * (a - r) + z * z;
    let (k, e) = ellip_agm(4.0 * a * r / p2);
    let c = MU_0 * i / (2.0 * PI * p2.sqrt());
    let bz = c * (k + (a * a - r * r - z * z) / q2 * e);
    let br = c * z / r * (-k + (a * a + r * r + z * z) / q2 * e);
    (br, bz)
}

/// Mutual inductance of coaxial loops (radii a, b, axial gap d) from Neumann's formula,
/// M = mu0 a b * int_0^pi cos(phi) / sqrt(a^2 + b^2 + d^2 - 2ab cos(phi)) dphi,
/// by the trapezoid rule over the full period.
pub fn neumann_coaxial(a: f64, b: f64, d: f64, n: usize) -> f64 {
    let dphi = 2.0 * PI / n as f64;
    let sum: f64 = (0..n)
        .map(|k| {
            let c = (k as f64 * dphi).cos();
            c / (a * a + b * b + d * d - 2.0 * a * b * c).sqrt()
        })
        .sum();
    0.5 * MU_0 * a * b * sum * dphi
}

/// Maxwell's closed form for coaxial loops, M = mu0 sqrt(ab) [(2/k - k) K - (2/k) E],
/// k^2 = 4ab / ((a + b)^2 + d^2), with [ellip_agm].
pub fn maxwell_mutual(a: f64, b: f64, d: f64) -> f64 {
    let p2 = (a + b) * (a + b) + d * d;
    let m = 4.0 * a * b / p2;
    let k = m.sqrt();
    let (kk, ee) = ellip_agm_c(m, ((a - b) * (a - b) + d * d) / p2);
    MU_0 * (a * b).sqrt() * ((2.0 / k - k) * kk - 2.0 / k * ee)
}

/// Kirchhoff's thin-ring self-inductance, uniform current in a round wire of radius b on a ring
/// of radius a, b << a: L = mu0 a (ln(8a/b) - 7/4).
pub fn kirchhoff_ring(a: f64, b: f64) -> f64 {
    MU_0 * a * ((8.0 * a / b).ln() - 1.75)
}

/// Wheeler (1928) multilayer air coil, L[uH] = 0.8 a^2 N^2 / (6a + 9b + 10c) with lengths in
/// inches: a mean radius, b length, c winding depth. Quoted to about 1% when the three
/// denominator terms are comparable.
pub fn wheeler_multilayer(a: f64, b: f64, c: f64, n: f64) -> f64 {
    let inch = 0.0254;
    let (a, b, c) = (a / inch, b / inch, c / inch);
    0.8 * a * a * n * n / (6.0 * a + 9.0 * b + 10.0 * c) * 1e-6
}

/// Wheeler (1928) single-layer solenoid, L[uH] = a^2 N^2 / (9a + 10l), inches. Quoted to about
/// 1% for l > 0.8a. It approximates the current sheet (Nagaoka), not round wire.
pub fn wheeler_single(a: f64, l: f64, n: f64) -> f64 {
    let inch = 0.0254;
    let (a, l) = (a / inch, l / inch);
    a * a * n * n / (9.0 * a + 10.0 * l) * 1e-6
}

/// Inductance of a uniform cylindrical current sheet (radius a, length l, N turns), the quantity
/// Nagaoka's coefficient tabulates: L = (N/l)^2 int_{-l}^{l} (l - |u|) M(a, a, u) du, integrated
/// numerically with [maxwell_mutual]. The log singularity at u = 0 is graded away with u = l s^2.
pub fn current_sheet(a: f64, l: f64, n: f64) -> f64 {
    let steps = 200_000;
    let ds = 1.0 / steps as f64;
    let mut sum = 0.0;
    for i in 0..steps {
        let s = (i as f64 + 0.5) * ds;
        let u = l * s * s;
        sum += (l - u) * maxwell_mutual(a, a, u) * 2.0 * l * s * ds;
    }
    // Both signs of u.
    (n / l).powi(2) * 2.0 * sum
}

/// Centre field of a thick solenoid with uniform current density (Montgomery 1969; the Fabry
/// form): B0 = mu0 J a1 beta ln[(alpha + sqrt(alpha^2 + beta^2)) / (1 + sqrt(1 + beta^2))],
/// alpha = a2/a1, beta = half-length / a1, J = N I / ((a2 - a1) * length).
pub fn thick_solenoid_centre(a1: f64, a2: f64, length: f64, ni: f64) -> f64 {
    let alpha = a2 / a1;
    let beta = 0.5 * length / a1;
    let j = ni / ((a2 - a1) * length);
    MU_0 * j
        * a1
        * beta
        * ((alpha + (alpha * alpha + beta * beta).sqrt()) / (1.0 + (1.0 + beta * beta).sqrt())).ln()
}
