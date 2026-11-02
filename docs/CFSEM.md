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
