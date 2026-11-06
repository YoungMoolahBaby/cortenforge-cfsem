# Toy stellarator (vacuum field, tabletop scale)

Modular stellarator coils as Fourier-series curves, their field from cfsem's straight-filament
Biot-Savart, field lines and Poincare sections to see whether they make nested magnetic surfaces,
and the coils swept into solids with cf-design for printing and winding. The target is a small
vacuum-field device (no plasma): 3D-printed coil formers wound with copper, run at low field,
checked with a Hall probe.

- cfsem: `flux_density_linear_filament` for the coils (each coil a closed polygon of straight
  filaments); later `inductance_linear_filaments_matrix_par` and the body-force functions.
- CortenForge: `cf_design::Solid::pipe_spline` sweeps each coil, `mesh_adaptive_par` meshes it,
  mesh-repair validates and mesh-io writes STL; later mesh-sdf for clearances, mesh-printability
  for the formers, sim-opt / sim-rl for the coil search.
- Ours: the Fourier coils (`src/coils.rs`), an RK4 field-line tracer in the toroidal angle, the
  period map's fixed points and the rotational transform (`src/trace.rs`).

## Modes

`cargo run --release -p stellarator -- <mode>`

- `check` (default): the coil model, cfsem's field for these coils and the tracer against known
  results, 21 checks (about 25 s). Exits 1 on any failure.
- `scan`: axis and rotational transform across coil shapes (field periods, ellipse, helical
  shift); how the baseline was picked (about 70 s).
- `trace`: the baseline's fixed points, 24 field lines for 600 periods, the iota profile and
  Poincare sections at phi = 0 and the half period; writes `out/poincare.svg` (about 80 s).
- `mesh`: the coils swept, meshed, snapped and checked; writes `out/mesh/coils.stl` (mm, 610k
  faces) and a render `out/mesh/coils.bmp` (about 3.5 min).
- `mesh-probe`, `sdf-probe`: the mesher findings below, reproduced on one coil.

## The baseline (`src/design.rs`)

A rotating ellipse with a helical shift: 3 field periods, 18 planar elliptical coils (6 a period)
in poloidal planes, R0 0.25 m, semi-axes 100 x 40 mm, the long axis turning by pi each period
(alpha = n_fp phi / 2), coil centres on a helix (R, Z) = (0.25 + 0.02 cos 3phi, -0.02 sin 3phi) m.
Each coil is an order-1 Fourier curve (the optimizer will add higher orders), cut into 128
straight filaments. 6.94 kA-turns a coil gives a toroidally averaged 0.1 T at R0. Coils sit
half a spacing off phi = 0, so phi = 0 and the half-period plane are planes of stellarator symmetry.

## What `check` found (2026-11-06)

- Fourier coils: an order-1 coil's length (trapezoid rule on |x'(t)|) equals the ellipse
  perimeter (Gauss-Kummer series) to rounding; sampled points lie on the ellipse in its plane to
  1.3e-15; the analytic tangent matches a central difference to 2e-11.
- cfsem's filaments against its own circular filament (a circular coil from the same code): max
  relative field error 1.82e-3 / 1.14e-4 / 7.10e-6 at 64 / 256 / 1024 segments, order 2.000.
- Ampere: a toroidal loop gives mu0 N I to 3.3e-13 and 3.6e-13 (two loops inside every coil); a
  poloidal loop between coils links nothing (1.8e-17 of B0 2 pi rho).
- div B and curl B by central differences, relative to |B| / 5 cm: 6e-9 and 2.8e-7, the same
  with and without cfsem's 4 mm wire radius. (At a 1e-4 m step curl B read 2.8e-5: truncation
  error 3 cm from a coil, falling as h^2.)
- Field-period symmetry 7.4e-15, stellarator symmetry (B_R odd, B_phi and B_Z even under
  (phi, Z) -> (-phi, -Z)) 3.6e-15.
- RK4 over one period: 2.6e-6 / 7.6e-8 / 2.0e-9 m at 60 / 120 / 240 steps against 1920 (faster
  than order 4 over a whole period); 120 steps a period is used.
- Circular planar coils (no ellipse, no shift): iota -2.3e-13, as it must be.
- The baseline axis crosses phi = 0 at Z = -4.6e-8 m (stellarator symmetry, to the RK4 error).
  det J of the period map at the axis is 1 to 2.4e-6 at 120 steps and 1.5e-8 at 240 (the map
  preserves the flux B_phi dR dZ; RK4 does so only to its step error).
- iota near the axis two ways: a line 2 mm off the axis over 600 periods (weighted Birkhoff
  average) 0.104991, the period map's Jacobian (n_fp acos(tr/2) / 2 pi) 0.105027, 3.4e-4 apart.

## Picking the coils (`scan`, 2026-11-06)

Planar coils only make a little rotational transform. At 18 coils:

| field periods | ellipse (mm) | helix shift R, Z (mm) | axis R (m) | iota on axis | iota 50 mm out |
|---|---|---|---|---|---|
| 3 | 100 x 50 | 0, 0 | 0.2613 | 0.0609 | 0.0476 |
| 3 | 100 x 50 | +20, -20 | 0.2663 | 0.0611 | 0.0509 |
| 3 | 100 x 40 | 0, 0 | 0.2606 | 0.0988 | lost |
| 3 | 100 x 40 | +20, -20 | 0.2644 | 0.1050 | 0.0755 |
| 3 | 100 x 40 | -20, +20 | 0.2529 | 0.1044 | lost |
| 6 | 100 x 50 | +20, -20 | 0.2597 | 0.0491 | 0.0530 |

- The ellipse does the work: 100 x 50 mm gives iota 0.06, 100 x 40 mm 0.10. Six periods with 3
  coils each is no better than three with 6: the ellipse is too coarsely sampled.
- The helical shift adds little transform (0.0988 -> 0.1050) but sets how far out lines stay
  confined, and its handedness matters: (+20, -20) keeps lines to 50 mm, (-20, +20) loses them.
- The period map has more than one fixed point (up to 3 found from 7 starting guesses). The first
  scan picked the one with the smallest Newton residual and once landed on the wrong one (R 0.185,
  "iota" -0.07). The axis is now the elliptic fixed point (|trace J| < 2) nearest the coil centre.

## Field lines (`trace`, 2026-11-06)

Fixed points at phi = 0: the axis (R 0.26440 m, elliptic, trace 1.9518), a hyperbolic point at
(0.22966, -0.0492) m, and an elliptic point at R 0.16781 m. The axis is at R 0.25069 m on the
half-period plane. 24 lines from the axis outward along the phi = 0 midplane, 600 periods (200
toroidal turns) each, 74 s:

| start (mm) | iota | height/width phi = 0 | half period |
|---|---|---|---|
| 2.5 | 0.10497 | 0.704 | 1.189 |
| 10 | 0.10422 | 0.705 | 1.179 |
| 20 | 0.10171 | 0.705 | 1.144 |
| 30 | 0.09718 | 0.702 | 1.081 |
| 40 | 0.09013 | 0.690 | 0.968 |
| 50 | 0.08266 | 0.652 | 0.760 |
| 52.5 | 0.08224 | 0.638 | 0.716 |
| 55 | lost in the first period | | |
| 57.5 | 0.08231 | 0.603 | 0.765 |
| 60 | lost after 12 periods | | |

- Nested surfaces out to about 52 mm from the axis on the outboard midplane; iota falls from
  0.105 to 0.082 (some magnetic shear). Birkhoff and plain averages agree to 1e-4 or better.
- The surfaces turn with the coils: wide at phi = 0 (height/width 0.70), tall at the half period
  (1.19) near the axis. Further out the outboard side draws out into a point (a teardrop), and in
  the plot that hides the rotation.
- The edge is ragged where iota crosses 1/12 (0.0839 -> 0.0822 between 45 and 52.5 mm): the
  12-island resonance. The line at 57.5 mm stays confined but outside it, past the separatrix.
- Cost: 6.9M field evaluations (2,304 filaments, one point each) in 74 s on 12 threads, about
  0.13 ms an evaluation a thread.

## Coils to STL (`mesh`, 2026-11-06)

Each coil: 128 points on its Fourier curve into `Solid::pipe_spline` (4 mm radius, in mm),
`mesh_adaptive_par` at 1 mm cells, then two Newton steps per vertex onto the solid's zero level
set along its SDF gradient. Checked against the true tube (distance from the Fourier centreline
minus r, 8192 samples) and pi r^2 L.

| coil 1 (460.3 mm) | faces | surface error (mm) | volume vs pi r^2 L | time |
|---|---|---|---|---|
| uniform `mesh`, 32 points, 1 mm | 34,894 | | -1.22e-2 | 305 s |
| octree, 64 points, 1 mm | 34,892 | -1.025 to +0.860 | -1.37e-3 | 3.9 s |
| octree, 128 points, 1 mm | 34,888 | -1.025 to +0.860 | -1.41e-3 | 8.0 s |
| octree, 128 points, 1 mm, snapped | 34,888 | -0.000 to +0.012 | -8.71e-3 | 8.5 s |
| octree, 128 points, 0.5 mm, snapped | 139,688 | -0.000 to +0.015 | -2.04e-3 | 35.3 s |

- `pipe_spline` is open-curve only; the control points are wrapped so the closed coil has no
  kink (see `docs/CORTENFORGE.md`). Its SDF matches the true tube to 0.06 mm at 64 points.
- The octree mesher is about 180x faster than the uniform one (1.7 s against 305 s at 32 points),
  but places vertices up to 0.9 of
  a cell off the surface (median 0.0075 mm at 0.5 mm cells). Snapping fixes that; what is left in
  the volume is the flat facets cutting chords across the round tube (-0.87% at 1 mm, -0.20% at 0.5 mm).
- No weld: the octree output is clean, and `weld_vertices` made it non-manifold.
- All 18 coils: watertight and manifold, surface within 0.015 mm, volume at worst -1.11% (the
  faceting), 143 s. Closest neighbouring centrelines 62.2 mm: 54.2 mm between conductors.

## Plan

1. `check` (done): coils, field, tracer.
2. Toy coils to surfaces (done): scan, axis, iota profile, Poincare sections.
3. Coils to STL with cf-design (done): swept, meshed, snapped, checked.
4. Next: score the surfaces (volume inside the last closed surface, iota away from low-order
   rationals, a smooth edge), penalize curvature, coil-to-coil distance and forces (cfsem body
   force), and let sim-opt / sim-rl move the Fourier coefficients; then printable coil formers,
   inductance and a Hall-probe field map to check a build.
