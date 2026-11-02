//! Warm-up 1: coil designer. For now, a smoke test that cfsem works:
//! the field at the centre of a Helmholtz pair against the textbook formula.

use cfsem::MU_0;
use cfsem::physics::circular_filament::flux_density_circular_filament_scalar;

fn main() {
    let (radius, turns, current) = (0.10_f64, 100.0, 1.0); // m, turns, A
    let ni = turns * current; // A-turns
    // Helmholtz: two loops, spacing = radius, centred on z = 0.
    let loops = [(radius, -radius / 2.0, ni), (radius, radius / 2.0, ni)];
    let bz_at = |r: f64, z: f64| -> f64 {
        loops.iter().map(|&l| flux_density_circular_filament_scalar(l, (r, z)).1).sum()
    };

    let b0 = bz_at(0.0, 0.0);
    let textbook = (4.0_f64 / 5.0).powf(1.5) * MU_0 * ni / radius;
    println!("Helmholtz pair, R = {radius} m, NI = {ni} A-turns");
    println!("  Bz(centre) cfsem    = {:.6e} T", b0);
    println!("  Bz(centre) textbook = {:.6e} T", textbook);
    println!("  relative error      = {:.2e}", (b0 - textbook).abs() / textbook);
    for frac in [0.1, 0.2, 0.3] {
        let dz = frac * radius;
        println!(
            "  uniformity at z = {:.2} R: {:+.3e}, at r = {:.2} R: {:+.3e}",
            frac, bz_at(0.0, dz) / b0 - 1.0, frac, bz_at(frac * radius, 0.0) / b0 - 1.0
        );
    }
}
