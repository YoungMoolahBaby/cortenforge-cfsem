# Coil designer (warm-up 1)

Optimize Helmholtz or Maxwell coils for a uniform field (Maxwell coils: a uniform gradient).
The gentlest way to learn cfsem.

- cfsem: circular filaments, `flux_density_circular_filament_scalar` and friends.
- CortenForge: optimize radius, spacing and turns (sim-opt), export coil formers to print.

## Modes

`cargo run --release -p coil-designer -- <mode>`

- `check` (default): cfsem and the coil model against known results, 58 checks. Exits 1 on any
  failure.
  - 1-8, cfsem's circular filaments: on-axis loop, off-axis field vs Biot-Savart quadrature and an
    exact (AGM) closed form, dipole far field, div B and curl B, Helmholtz and Maxwell-pair
    textbook values, mutual inductance vs Neumann, thin-ring self-inductance vs Kirchhoff.
  - 9-12b, the coil model (`src/coil.rs`): pack inductance vs independent kernels, Wheeler's
    multilayer and single-layer formulas, an exact current-sheet integral, Montgomery's
    thick-solenoid centre field, and the uniform-sphere search vs the Helmholtz series.
  - 13, speed.
- `pack`: a tabletop Helmholtz pair wound from AWG 22, its electrical numbers, and its uniform
  volume against the ideal thin-loop pair, with a spacing scan.

Reference formulas written independently of cfsem live in `src/refs.rs`.

## The coil model

Each turn of a layer-wound rectangular pack is one cfsem filament at its centre (square packing,
pitch = wire diameter over enamel, helix ignored). Field: ideal filaments, exact outside the wire.
Inductance: mutual terms from the ideal-filament potential (exact between round wires), self terms
from cfsem's finite-wire model averaged over the wire section. Resistance from the IACS copper
resistivity at 20 C. AWG diameters by definition; the 0.06 mm enamel build is an assumption.

## What it found (2026-11-02)

`check`:
- All 58 pass. cfsem's off-axis field is 0.8-1.3e-8 relative off the exact closed form; that's
  its elliptic-integral curve fits, not the physics (the exact form and the quadrature agree to
  1e-14). Mutual inductance agrees with Neumann to 1.6e-9.
- Helmholtz flatness follows the textbook series: B(z)/B0 - 1 = -1.842e-7 at z = 0.02 R
  (series -1.843e-7). Maxwell-pair gradient 8.0587e-3 T/m for NI = 100 A, R = 0.1 m, matching
  0.6413 mu0 NI/R^2 to 1e-5 (the 4-digit constant's rounding).
- Coil model: inductance within 1.4e-6 of independent kernels and 0.54% of Wheeler's multilayer
  formula (quoted to ~1%); thick-solenoid centre field within 2e-5 of Montgomery. A single-layer
  round-wire coil comes out 1.0% below the ideal current sheet: the sign expected for round wire
  (Rosa's correction), and the order of mu0 a N / L = 1.3%.
- Speed: 1M points x 2 loops in 93-105 ms serial, 23-31 ms with `_par`.

`pack`, two coils of 100 turns AWG 22 (10 x 10), mean radius 50 mm, pack 7.0 x 7.0 mm:
- Series: 62.8 m of wire, 183 g copper, R = 3.33 ohm, L = 4.08 mH, L/R = 1.23 ms.
- B0 = 1.798 mT/A, 3.2e-4 below the thin-loop formula at the mean radius. 1 mT takes 0.556 A,
  1.85 V, 1.03 W.
- Uniform sphere (|B - B0|/B0 <= eps), thin loops vs real packs at spacing = mean radius:

  | eps | thin loops | real packs |
  |---|---|---|
  | 1e-2 | 15.69 mm | 15.64 mm |
  | 1e-3 | 8.65 mm | 8.55 mm |
  | 1e-4 | 4.84 mm | 4.66 mm |
  | 1e-5 | 2.72 mm | 2.40 mm |

- The finite section costs uniformity, mostly at tight tolerances. Spacing can buy it back: at
  1.005 R the real pair's 1e-4 sphere grows to 5.31 mm, past the thin pair at d = R (a small
  second-order term traded against the fourth-order one). The thin pair hasn't been given the
  same chance yet; step 2 optimizes both.

## Next

1. Optimize spacing and pack shape for the largest uniform sphere with sim-opt, checked against
   a brute-force scan, comparing thin and real pairs fairly (both optimized).
2. A printable coil former with cf-design, exported and checked for printability.
