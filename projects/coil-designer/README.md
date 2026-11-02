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
- `design`: the least-power pair meeting a robust uniformity spec (about 2 minutes), then its
  former.
- `former [turns per layer] [spacing mm] [cell mm] [mc|dc|adaptive]`: the printed coil former
  alone (defaults: the chosen design, 0.25 mm marching cubes). Writes `out/former_half.stl` and
  `out/former_washer.stl` (git-ignored).

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

## `design` (2026-11-02)

The power-aware objective: least power at 1 mT, subject to a 10 mm diameter sphere uniform to
1e-4 that holds for any spacing error within +-0.1 mm (an FDM tolerance), packs at least 5 mm
apart. Each coil is exactly 100 turns of AWG 22 at 50 mm mean radius (a partial last layer, so
shapes compare at equal copper). For each axial turn count (4-100) it finds the spacing window
that meets the spec and centres the design in it. 97 shapes, 103 s.

- 96 of 97 shapes have a window of at least 0.2 mm (4 turns/layer misses by 0.01 mm). Windows
  are tight for compact packs (0.21-0.31 mm for 5-24 turns/layer) and widen for long ones.
- Power is nearly flat across compact packs (1.010-1.040 W for 5-16 turns/layer) and climbs for
  long ones (1.225 W at 2 layers, 1.865 W at 1 layer); the long ones buy more margin, not more
  spec.
- Chosen: 7 turns/layer x 15 layers (last layer 2 turns), 4.93 x 10.56 mm, r_in 44.72 mm, at
  49.542 mm spacing (0.9908 R; window 49.43-49.65 mm). The sphere is 5.38 mm at the centre and
  >= 5.02 mm across +-0.1 mm. 1.811 mT/A: 1 mT at 0.552 A, 1.83 V, 1.010 W; R 3.311 ohm,
  L 3.975 mH, 62.5 m of wire.

## The former (2026-11-02)

First try, one piece printed upright: printability flagged three 12 mm annular channel ceilings
at 90 deg, which would need support inside the winding channels. Redesign, two part types:

- **Half** (print 2): tube (bore 84.44 mm, wall 2.5 mm), winding floor flange (113.56 mm) with a
  40 deg support cone under it, and a bed collar with three blind dowel holes for 1.75 mm
  filament, pointed so their ends need no bridge. A 1.2 mm radial slot through the floor flange
  and cone (at 60 deg, between dowels) lets the start lead, which begins on the innermost layer
  against the floor, drop out without crossing the winding. Printed standing on its mid-plane
  face; the halves butt there.
- **Washer** (print 2): flat ring that slides over the tube end to close the channel, with a
  1.2 mm notch at its rim for the finish lead; it turns freely, so the notch lines up with
  wherever the winding ends.
- The winding rests on the floor, so the spacing is 2 x (floor + pack/2). The floor is snapped to
  0.1 mm layers: 22.300 mm, giving 49.527 mm (-0.015 mm from the design, inside the window).
  At 0.2 mm layers a slicer would have put it at 22.4, +0.185 mm on the spacing, past the whole
  tolerance. The pack width assumes 0.704 mm wire over enamel; 7 turns of wire 0.01 mm thicker
  move each coil 0.035 mm, so measure the wire and set `ENAMEL_BUILD` before printing.
- Mesh: marching cubes at 0.25 mm, welded at 1 um; worst vertex 0.062 mm off the solid; SDF
  probes on the floor, winding surface, bore and mid-plane faces read 0.0000 mm. The washer is
  "ready for printing". The half's remaining flags are 0.0 mm^2 slivers on the support cone at
  mesh-grid planes (details in `docs/CORTENFORGE.md`); the orientation search keeps both parts as
  modelled.
- SDF probes confirm the slots: +0.60 mm (open, half the slot width) on their centre lines,
  solid one slot width aside, in the flange, the cone and the washer.
- Not yet modelled: a way to fix the washers (glue for now); the leads are not in the field
  model (twist them).

## Next

1. Print, wind and measure B0 and the field profile against the model (a hardware step).
2. A Maxwell gradient pair through the same pipeline.
