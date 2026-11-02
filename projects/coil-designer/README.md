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
- `optimize`: spacing and pack shape for the largest uniform sphere: a brute-force grid, then
  CortenForge's gradient-free optimizers on the same problem (about 6 minutes).

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
  second-order term traded against the fourth-order one). `pack` scans in 0.0025 R steps;
  `optimize` finds the peak at 1.0039 R, 5.517 mm, equal to the thin pair at its own best spacing.

## `optimize` (2026-11-02)

Largest 1e-4 sphere for a pair of ~100-turn AWG 22 coils at 50 mm mean radius, over spacing
(0.9-1.6 R) and layer count (1-25, turns per layer = round(100 / layers)), with at least 5 mm
between the packs (two 2.5 mm former walls, an assumption). Run 3 of 3; the numbers below are
from it.

Brute force (10,800 evaluations, 90 s), best spacing per layer count, a few rows:

| Pack (layers x turns) | Best spacing | Gap | 1e-4 sphere | Field | Power for 1 mT |
|---|---|---|---|---|---|
| 1 x 100 (0.7 x 70.4 mm) | 1.5729 R | 8.3 mm | 8.177 mm | 1.336 mT/A | 1.87 W |
| 2 x 50 | 1.1508 R | 22.3 mm | 6.210 mm | 1.649 mT/A | 1.22 W |
| 4 x 25 | 1.0393 R | 34.4 mm | 5.688 mm | 1.757 mT/A | 1.08 W |
| 10 x 10 (7.0 x 7.0 mm) | 1.0039 R | 43.2 mm | 5.517 mm | 1.794 mT/A | 1.03 W |
| 25 x 4 | 0.9719 R | 45.8 mm | 5.371 mm | 1.825 mT/A | 1.00 W |
| thin loops (reference) | 1.0030 R | | 5.517 mm | | |

- Uniformity alone picks a split solenoid: one layer, 70 mm wide packs 8 mm apart, a sphere 48%
  larger (3.3x the volume) than the square pack, for 1.8x the power per mT. Field per amp and
  power vary a few % row to row because rounding changes the turn count (92-110).
- With its spacing optimized, the square 10 x 10 pack matches ideal thin loops (5.517 mm both):
  the finite section costs nothing at 1e-4 once spacing is tuned.
- The peak is sharp: the 10 x 10 pair is 5.52 mm at 50.195 mm spacing and 4.66 mm at 50.00 mm
  (`pack`). A 0.2 mm build error costs 16%, so former tolerances matter more than the section.

CortenForge's optimizers (3 seeds each; shipped = `best_artifact`), median gap to the grid:

| Problem | SA | RicherSA | PT | CEM |
|---|---|---|---|---|
| spacing + layers, 400 evals | -30.9% | -34.3% | -33.3% | -33.1% |
| spacing + layers, 1600 evals | -28.8% | -34.3% | -31.2% | -32.8% |
| spacing only at 10 x 10, 200 evals | -0.23% | +0.14% | -0.50% | -4.1% |

- On one parameter they work, RicherSA best (its adaptive step beats the 5e-5 R grid). On two
  coupled parameters none found the single-layer design: each layer count has its own narrow
  best spacing, so local moves fall off the ridge. For design work here, a grid over the integer
  with a 1-D search inside is the tool.
- CEM ships its unevaluated elite mean, always below its best evaluated design (up to 7.7%).
- Details and file:line references in `docs/CORTENFORGE.md`.

Earlier runs: run 1 searched 4-25 layers and 0.95-1.10 R and found the optimum on the box edge
(4 layers), which prompted the wider box and the gap constraint; run 2 rounded the layer count,
giving the end values half-width bins, fixed in run 3 (no change in outcome).

## Next

1. A printable coil former with cf-design for the chosen design, exported and checked for
   printability, with the spacing tolerance in mind.
2. Optionally, a power- or field-constrained objective (uniform sphere per watt), since
   uniformity alone drifts to long, inefficient coils.
