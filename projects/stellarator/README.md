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
- `trace`: the baseline's fixed points, 36 field lines (2.5-90 mm out) for 600 periods, the iota
  profile and Poincare sections at phi = 0 and the half period; a line is lost within 10 mm of a
  coil centreline or 120 mm from the coil centre. Writes `out/poincare.svg` (about 2 min).
- `mesh`: the coils swept, meshed, snapped and checked; writes `out/mesh/coils.stl` (mm, 604k
  faces) and a render `out/mesh/coils.bmp` (about 2 min).
- `mesh-probe`, `sdf-probe`: the mesher findings below, reproduced on one coil.
- `opt-probe`: the baseline through the coil search's score, line by line.
- `optimize [generations] [population] [seed]`: the coil search (CEM, default 24 x 24, about
  1.5 h); writes `out/opt/evals.csv` and `out/opt/best.txt`.
- `best [design.txt] [--no-trace] [--no-mesh]`: a saved design (default `out/opt/best.txt`;
  run 1's winner is kept in `designs/run1.txt`) at full resolution: `trace` into
  `out/poincare_best.svg` and `mesh` into `out/mesh_best/` (about 4 min).

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
- Since the coil search, `trace` uses the search's loss rule (within 10 mm of a coil centreline,
  or 120 mm from the coil centre) and runs lines to 90 mm. Rerun that way the baseline gives the
  same table: confined to 52.5 mm without a gap, 55 mm lost (it nears a coil), 57.5 mm outside
  the separatrix, beyond lost.
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
- All 18 coils: one field period (6 coils) meshed and turned into the other two, since the
  octree mesher gave 8-10 non-manifold edges on 4 of the 12 rotated copies of the optimized coils
  when they were meshed in place (`docs/CORTENFORGE.md`). Watertight and manifold, surface within
  0.015 mm, volume at worst -1.11% (the faceting), 57 s. Closest centrelines 62.2 mm: 54.2 mm
  between conductors.

## Coil search (`optimize`, 2026-11-06)

**Design space** (`src/search.rs`). The three coils of a half period are free; the other fifteen
follow by stellarator symmetry (rotate by pi about x and run backwards) and the field period, so
every candidate keeps both symmetries. Built from the baseline's three coils, the set's field
matches the baseline's to 3.3e-14. Each free coil is the baseline coil plus a displacement in
its own frame (e_R, e_phi, e_Z) to Fourier order 2: 45 numbers, each a fraction of 10 / 15 / 8 mm
for orders 0 / 1 / 2. The e_phi modes bend coils out of their planes.

**Score.** Twelve lines start 5-60 mm out from the axis and run 300 periods (64 filaments a coil,
120 RK4 steps a period). A line is good if it stays 10 mm from every coil centreline (4 mm
conductor plus former) and the Birkhoff iota of its two halves agree to 1e-3. Score = the
outermost good line in an unbroken run from the axis (mm) + 50 x the lowest iota on those lines
(capped at 0.25), minus penalties below a 12 mm bend radius or 30 mm between centrelines. The
baseline scores 49.26 (45 mm, iota 0.105 to 0.085), about 10 s an evaluation.

Calibrating the score (`opt-probe`) took three fixes, and a fourth came after the run:
- 60 RK4 steps a period (fine for one period, 7.6e-8 m) let outer lines drift over 300 periods:
  their iota climbed from one half to the other (0.0919 -> 0.0941 at 40 mm) and good surfaces read
  as chaos. 120 steps matches the full-resolution trace; 64 or 128 filaments made no difference.
- 150 periods is too short: a half then covers about two poloidal turns and real surfaces at
  40-45 mm fail the halves test (1.0e-3 and 1.5e-3); at 300 periods they read 1.7e-4 and 7.8e-4,
  and the chaotic band at 50 mm 2.6e-3.
- 15 mm clearance cut the baseline's good region to 30 mm: surfaces do pass within 15 mm of the
  coils near the half period. 10 mm leaves room for a 4 mm conductor and a thin former.
- After the run: a line escaping outward never nears a coil, so the loss test now also drops lines
  120 mm from the coil centre. The best design is unaffected (its lines stay within 74 mm).

**Run 1.** sim-rl's CEM, 24 generations x 24 candidates, seed 1, noise 0.2 decaying 0.93, through
the dummy-MJCF wrapper from coil-designer; every evaluation logged (CEM's `best_artifact` is an
elite mean it never evaluated). 577 evaluations in 94 min; 28 (5%) had no elliptic axis. The
population mean rose from -6.0 to 72.1 and the best, 72.50, came at evaluation 262 (generation
11). That is the score's ceiling: the good region reached the last line (60 mm) and iota passed
the 0.25 cap, so run 1 ended limited by the score, not the coils.

**The winner at full resolution** (`best`; 128 filaments, 600 periods, lines to 90 mm, same loss
rule):

| | baseline | optimized |
|---|---|---|
| axis R at phi = 0 | 0.26440 m | 0.25231 m |
| iota on axis (Jacobian / traced 2.5 mm) | 0.10503 / 0.10497 | 0.26918 / 0.26921 |
| iota at the edge of the good region | 0.0822 (52.5 mm, falling) | 0.2927 (65 mm, rising) |
| confined without a gap, outboard midplane | 52.5 mm | 65 mm |
| surface height / width near the axis, phi = 0 / half | 0.70 / 1.19 | 0.68 / 1.88 |
| coil lengths (mm) | 460.3 x 3 | 470.5 / 451.7 / 486.8 |
| out of plane | 0 | up to 9.8 / 14.5 / 5.6 mm |
| tightest bend, closest centrelines | 16.0 mm, 62.2 mm | 12.1 mm, 52.0 mm |

- Iota is 2.6x the baseline's and now rises outward (reversed shear). The good region is 25%
  wider; past 65 mm lines are lost within 3 periods by nearing a coil.
- Dashed rings in `out/poincare_best.svg` are island chains where iota passes simple ratios
  (likely 3/11 = 0.2727 near 27.5 mm and 2/7 = 0.2857 near 55-57.5 mm).
- The bend limit is active: the tightest bend sits at 12.1 mm against the 12 mm floor.
- Coils mesh watertight and manifold (one period meshed and turned), surface within 0.015 mm,
  `out/mesh_best/coils.stl` (625k faces).

## Plan

1. `check` (done): coils, field, tracer.
2. Toy coils to surfaces (done): scan, axis, iota profile, Poincare sections.
3. Coils to STL with cf-design (done): swept, meshed, snapped, checked.
4. Coil search (run 1 done): CEM on 45 Fourier coefficients, iota 0.105 -> 0.269, good region
   52.5 -> 65 mm.
5. Next: a score that can't saturate (lines out to 90 mm, iota cap raised or replaced by a
   target away from low-order rationals, island width), then more runs and seeds and sim-opt's
   annealers for comparison; coil forces (cfsem body force) and inductance; printable coil
   formers; a Hall-probe field map to check a build.
