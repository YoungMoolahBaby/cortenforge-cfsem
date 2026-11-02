# cfsem 14.0.0: end-user findings

Friction found while using the published crate as an end user. CortenForge's findings live in
`CORTENFORGE.md`. Log it while building, each entry starting with its type (bug, docs, API, perf,
setup, works) and the date, with file:line in the crate source
(`~/.cargo/registry/src/*/cfsem-14.0.0/`).

## cfsem 14.0.0

- works, 2026-11-02: `physics::circular_filament::flux_density_circular_filament_scalar` gives the
  Helmholtz centre field to 3e-12 of the textbook formula (coil-designer smoke test).
- works, 2026-11-02: the circular-filament field passes every textbook check in
  `cargo run --release -p coil-designer -- check` (43 of 43): on-axis loop, off-axis field against
  Biot-Savart quadrature and an exact (AGM) closed form, dipole far field, div B and curl B by
  finite differences, Helmholtz centre field and its -(144/125)(z/R)^4 flatness, Maxwell-pair
  gradient 0.6413 mu0 NI/R^2. Mutual inductance from `flux_circular_filament` with a 1 A source
  matches Neumann's formula to 1.6e-9.
- works, 2026-11-02: the finite-wire near-field model (Hurwitz et al.) gives thin-ring
  self-inductance: the cross-section mean of `flux_circular_filament_scalar` matches Kirchhoff's
  mu0 a (ln(8a/b) - 7/4) to 8e-6 at b/a = 0.01 and 7e-5 at b/a = 0.03.
- docs, 2026-11-02: the field functions say "Near-exact formula (except numerically-evaluated
  elliptic integrals)" (physics/circular_filament.rs:454, 1220) but not what that costs. The
  elliptic integrals are curve fits (math.rs:113, 142: 2e-8 absolute; math.rs:172: 4e-9), and the
  off-axis field comes out 0.8-1.3e-8 relative off an exact closed form. On the axis
  `flux_density_circular_filament_scalar` is off by a constant 3.2e-12 (the E fit at m = 0), while
  `flux_density_circular_filament_on_axis` is exact. Fine for coil design; worth one line on the
  field functions for anyone mapping ppb-level homogeneity.
- docs, 2026-11-02: there is no function named for mutual inductance between coaxial loops. It's
  `flux_circular_filament` with a 1 A source, observed at the second loop; the docs say "Flux is
  2 pi R A_phi" (physics/circular_filament.rs:61, 95) but never "this is M per amp".
- API, 2026-11-02: the Rust crate lacks the coil-level helpers the Python package has: Lyle's and
  Wien's self-inductance, `self_inductance_axisymmetric_coil`, `mutual_inductance_of_cylindrical_coils`
  (cfsem/__init__.py:248, 463, 562, 620, written in Python). A Rust user filamentizes winding packs
  and sums self-inductance themselves: about 60 lines in coil-designer (src/coil.rs `inductance`,
  `ring_self_inductance`).
- API, 2026-11-02: `flux_circular_filament_scalar` requires a nonzero wire radius
  (physics/circular_filament.rs:115) while the vector form takes zero for an ideal filament, so the
  scalar ideal-loop flux means scaling `vector_potential_circular_filament_scalar` by 2 pi R by hand.
- perf, 2026-11-02: 50 ns per loop-point serial (1M points x 2 loops in 101 ms);
  `flux_density_circular_filament_par` takes 31 ms on 6 cores / 12 threads, a 3.3x speedup.
- works, 2026-11-02: a turn-by-turn coil model on cfsem's kernels (one filament per turn, mutual
  terms from the ideal-filament potential, self terms from the finite-wire model) checks out:
  inductance matches independent kernels (Kirchhoff + Maxwell's closed form) to 1.4e-6 and
  Wheeler's multilayer formula to 0.54%; the centre field of thick solenoids matches Montgomery's
  uniform-current-density formula to 1.4-1.9e-5 (`coil-designer -- check`, checks 9-12).
- works, 2026-11-02: inside the wire, the finite-thickness potential is a cubic in the local
  coordinates (physics/circular_filament.rs, `t <= 1` branch of
  `vector_potential_circular_filament_finite_thickness_scalar_near`), so a 4 x 16 Gauss grid
  averages it exactly. Reading the source made that quadrature choice easy; the docstring could
  say it.
- works, 2026-11-02: the Grad-Shafranov operators check out (fusion-visualizer `check`, 16 of
  16): `gs_operator_order2` / `_order4` converge at 2.00 / 4.00 on a smooth test function, and
  both reproduce degree-4 polynomials, including a Solov'ev equilibrium solved as a fixed-boundary
  problem, to rounding (the 2nd-order one is in flux form, so it is exact there too). A coil's
  flux solved on a 161 x 161 grid matches `flux_circular_filament` to 4e-6. B from Psi converges
  to `flux_density_circular_filament` at order 4.17, and Delta* of cfsem's vacuum flux is zero to
  3.9e-8: flux and field conventions agree.
- API, 2026-11-02: the Rust crate stops at the operator. The solver (`flux_solver`, a scipy
  sparse LU), the 4th-order gradient and B-from-Psi live only in Python (cfsem/flux_solver.py:40,
  107, 144). A Rust user brings a sparse LU (faer here, already in the tree via CortenForge) and
  writes the gradient: about 150 lines in fusion-visualizer src/gs.rs.
- docs, 2026-11-02: the Rust operator docs (physics/gradshafranov.rs:44-61, 113-147) don't say
  how to use it: that Psi is total flux (2 pi R A_phi, Wb), so the right side is
  -2 pi mu0 R J_phi; that node (i, j) is row i * nz + j; or that boundary rows are identity rows
  taking the boundary Psi. The first is in the Python docstring (flux_solver.py:149), the other
  two only in the source.
- works, 2026-11-02: `flux_circular_filament_par` makes a good free-boundary Green's function.
  A Picard-iterated free-boundary equilibrium (fusion-visualizer `shafranov`, `shape`) takes the
  grid edge flux of the plasma from it each iteration; at convergence the grid solution matches
  the Green's sum for the same current to 3.4e-7 at 12,205 vacuum nodes. The vertical field the
  solver needs to centre a circular plasma matches Shafranov's formula to -2.71 / -1.38 / -0.22%
  at a/R0 = 0.3 / 0.2 / 0.1, shrinking with a/R0 as the large-aspect-ratio formula should.
- works, 2026-11-02: `flux_density_circular_filament_cartesian_scalar` (a loop at any position and
  orientation) agrees with an m-segment polygon of `flux_density_linear_filament_scalar` at
  1.7e-3 / 1.0e-4 / 6.6e-6 for m = 64 / 256 / 1024, the expected 1/m^2: two independent kernels,
  one answer. Sixteen of them as TF coils satisfy Ampere's law (toroidal average of B_phi =
  mu0 N I / 2 pi R) to 1e-8 or better and div B = 0 to 2e-10, and trace 12 field lines x 100 turns
  (1.8M field evaluations, half of them through all 16 coils) in 0.6 s on 12 threads
  (fusion-visualizer `trace`).
- API (scope), 2026-11-02: no magnetic materials ("quasi-steady electromagnetics including
  filamentized approximations, Biot-Savart, and Grad-Shafranov", README.md:3). Iron-core
  electromagnets, the usual EMS maglev magnet, and steel keepers or shields can't be modelled; a
  flat, wide, unsaturated iron surface can, by images (maglev-train). That steered the maglev
  design to a permanent-magnet hybrid: NdFeB (mu_r ~ 1.05) is close to exact as equivalent surface
  currents, while an air-core coil at tabletop scale makes 0.049 N at 1 A and 5 mm (a 300 g cart
  would burn ~170 W).
- works, 2026-11-02: a disc magnet as equivalent current loops plus its iron image agrees with two
  `flux_density_dipole_scalar` dipoles (point_source/dipole.rs:20) to 1.5e-5 - 1.1e-4 at 0.14 to
  1 m; force between coaxial loops by Lorentz (B_r), by energy (flux) and by Maxwell's closed form
  agree to 1e-8 (maglev-train `check`). The batch oriented-loop field
  (`flux_density_circular_filament_cartesian_par`) maps a 4-magnet cart's 3D stray field
  directly.
- bug, 2026-11-02: `mutual_inductance_circular_to_linear_scalar` (physics/circular_filament.rs:1434)
  takes the segment midpoint wrongly. Its radius is `path_r + path_dr / 2.0` (:1453): the start
  radius plus half the segment's whole horizontal length, not half its change in radius, so only
  a segment pointing straight outward gets the right radius. Its angle is the mean of the two
  endpoint `atan2` values (:1455), which wraps wrongly across the -x axis. Measured against the
  same formula evaluated at the true (Cartesian) midpoint, for a loop of radius 10 mm
  (maglev-train `motor`, probe 0): a 10 mm tangential segment at r 20 mm is 40% low; 10 mm
  diagonal segments are 14% low going outward and 46% low going inward; a 0.2 mm segment
  crossing the -x axis comes out with the wrong sign. In use (coenergy thrust of a flat stator
  coil on magnet loops), the error converges only to first order with segment length (9% at
  1 mm, 2% at 0.25 mm) where the midpoint rule should be second order. Fix: take the midpoint in
  Cartesian coordinates, r = hypot(xm, ym), phi = atan2(ym, xm). With that one change the same
  thrust matches an independent Lorentz-force calculation to 4.0e-4 at 1 mm and 2.6e-5 at
  0.25 mm (second order).
- works, 2026-11-02: `flux_density_linear_filament` for flat stator coils (exact straight
  segments, mirrored in a steel rail) drives a linear-motor model whose Lorentz-force thrust on
  magnet loops agrees with the coenergy route to 2.6e-5 once the midpoint above is corrected
  (maglev-train `motor`).
- works, 2026-11-02: a second, guidance coil layer for the maglev cart (coils spanning both magnet
  rows) gives sideways force by Lorentz force that agrees with the coenergy route
  (true-midpoint flux) to 7e-5 at 0.25 mm pieces, and side effects of 1e-16 where symmetry says
  zero (maglev-train `guide`).
- bug, 2026-11-02: linear-filament self-inductance counts the wire's internal inductance twice
  over. `inductance_piecewise_linear_filaments` (physics/linear_filament.rs:74) integrates A.dl
  at points on the wire axis (:22-48), where the in-wire term
  `blend = 1/2 (1 - frac^2)` (:1293) is at its peak: mu0/4pi per metre. The energy (and DC)
  value is the cross-section average, mu0/8pi per metre. The docs say "uniform current
  distribution; low frequency", so the result should be the DC one. A 50 mm ring of 1 mm wire
  (1024-gon) gives 325.76 nH against 310.06 nH from Wien's formula, which cfsem's own Python
  package ships (`self_inductance_circular_ring_wien`): +5.1%, and the excess (15.705 nH) is
  mu0/8pi x length (15.708 nH). A 100 x 50 mm rectangle: +5.6% against Grover's partial
  inductances. cfsem's own test (test/test_electromagnetics.py:247) compares a 13-turn coil
  with Lyle at 5% tolerance, where turn-to-turn terms swamp it. Workaround: subtract
  mu0/8pi x length; that matches Wien to 1.1e-5 and Grover to 7.6e-4. Mutual terms are
  unaffected (wireless-charging `check`).
- works, 2026-11-02: the linear-filament kernels otherwise check out (wireless-charging
  `check`, 9 of 9): mutual inductance of two offset squares against an independent
  Gauss-Legendre Neumann integral to 2e-10, both orders; polygons converge on
  `flux_circular_filament` at second order (1024-gons: 5.5e-6); `flux_density_linear_filament`
  matches the Biot-Savart finite-wire closed form to 1e-15; loop images satisfy the
  conductor and iron boundary conditions to 1e-16.
- works, 2026-11-02: `inductance_piecewise_linear_filaments` sizes an inductive power link. A
  pickup turn over a 20 m two-wire track matches the 2D flux to 4.5e-6, and a 12-turn spiral
  pickup's L2 moves 2e-5 between 5 mm and 2.5 mm pieces; a whole design sweep runs in about 8 s
  (wireless-charging `link`).
- API (scope), 2026-11-02: no eddy currents in conducting or magnetic plates, so a steel rail at
  85 kHz (skin depth 20-70 um) is only bounded by images (iron +1, conductor -1), and a ferrite
  pickup backing can't be modelled. Workaround that worked: a 2D spectral image weight
  R(k) = (mu_r k - gamma) / (mu_r k + gamma) for the long track, its ratio to the iron case
  scaling cfsem's 3D image terms (checked to 5e-8 against the closed-form images, and 0.7%
  against cfsem's 3D term). A complex image weight per term, or a conducting half-space Green's
  function, would cover this directly (wireless-charging `effects`).
- works, 2026-11-02: polygon loops through `inductance_piecewise_linear_filaments` handle offset
  and tilted coils: two circles 8 mm apart laterally, flat or tilted 10 deg, match an independent
  Neumann integral on exact circles to 5.6e-6 / 6.3e-6 at 1024 sides (second order). The Qi A11
  pad coil comes out 3.74 uH in free space and 6.99 uH over an ideal iron plane, against the
  spec's 6.3 uH +- 10% with real ferrite (wireless-charging `pad`).
