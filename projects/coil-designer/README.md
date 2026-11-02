# Coil designer (warm-up 1)

Optimize Helmholtz or Maxwell coils for a uniform field (Maxwell coils: a uniform gradient).
The gentlest way to learn cfsem.

- cfsem: circular filaments, `flux_density_circular_filament_scalar` and friends.
- CortenForge: optimize radius, spacing and turns (sim-opt), export coil formers to print.

## Modes

`cargo run --release -p coil-designer -- <mode>`

- `check` (default): cfsem's circular-filament API against known results, 43 checks: on-axis
  loop, off-axis field vs Biot-Savart quadrature and an exact closed form, dipole far field,
  div B and curl B, Helmholtz and Maxwell-pair textbook values, mutual inductance vs Neumann,
  thin-ring self-inductance vs Kirchhoff, and a speed line. Exits 1 on any failure.

## What `check` found (2026-11-02)

- All 43 pass. cfsem's off-axis field is 0.8-1.3e-8 relative off the exact closed form; that's
  its elliptic-integral curve fits, not the physics (the exact form and the quadrature agree to
  1e-14). Mutual inductance agrees with Neumann to 1.6e-9.
- Helmholtz flatness follows the textbook series: B(z)/B0 - 1 = -1.842e-7 at z = 0.02 R
  (series -1.843e-7). Maxwell-pair gradient 8.0587e-3 T/m for NI = 100 A, R = 0.1 m, matching
  0.6413 mu0 NI/R^2 to 1e-5 (the 4-digit constant's rounding).
- Speed: 1M points x 2 loops in 101 ms serial, 31 ms with `_par`.

## Next

1. Real winding packs: N turns in a rectangular section, filamentized; inductance (own Lyle or
   filament sums, since Rust cfsem has no coil-level helpers), resistance, power.
2. Optimize spacing and pack shape for the largest 1e-3 / 1e-4 uniform volume with sim-opt,
   checked against a brute-force scan.
3. A printable coil former with cf-design, exported and checked for printability.
