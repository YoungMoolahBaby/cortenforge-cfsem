# Fusion magnet visualizer (warm-up 2)

Trace field lines around a small tokamak or stellarator and show the plasma equilibrium shape.

- cfsem: circular filaments for the tokamak coils (TF, PF, central solenoid), linear filaments
  for stellarator coils, Grad-Shafranov (`cfsem::physics::gradshafranov`) for the equilibrium.
- CortenForge: coil and vessel geometry, mesh export of the coils and flux surfaces.
- Done when: field lines close on nested flux surfaces, and the equilibrium shape responds to the
  PF currents.

## Modes

`cargo run --release -p fusion-visualizer -- <mode>`

- `check` (default): cfsem's Grad-Shafranov operator and its flux and field conventions against
  known results, 16 checks. Exits 1 on any failure.
- `shafranov`: the free-boundary solver against Shafranov's vertical-field formula at
  a/R0 = 0.3, 0.2, 0.1 (about 5 s). `FV_DEBUG=1` prints every iteration.
- `shape`: a small tokamak's plasma shape across six PF coil settings; writes
  `out/shape_*.svg` (git-ignored).

## Grad-Shafranov plumbing (`src/gs.rs`)

cfsem's Rust crate ships the Delta* operator (2nd and 4th order) as sparse triplets; the solver
is ours: a regular (R, Z) grid, a faer sparse LU factorized once per grid, the right side
-2 pi mu0 R J_phi with boundary Psi from cfsem's filament flux, and B from Psi by 4th-order
differences. Psi is total poloidal flux (2 pi R A_phi, Wb), cfsem's convention throughout.

## What `check` found (2026-11-02)

- Delta* on a smooth function (R^2 cos R exp(-Z^2)) converges at 2.00 (order 2) and 4.00
  (order 4). Both stencils reproduce R^4, R^2 Z^2 and the Solov'ev equilibrium
  Psi = R^2 Z^2 + (a/4)(R^2 - R0^2)^2 to rounding (<= 7.7e-11), and the solved magnetic axis sits
  at (R0, 0) exactly.
- A coil block's flux solved on the grid matches cfsem's filament flux for the same source to
  2.7e-4 / 9.5e-5 (order 2, 81 / 161 nodes a side) and 1.3e-5 / 4.3e-6 (order 4). It converges at
  about order 1.5 either way: the source is point currents at nodes, not the stencil, sets it.
- B from Psi converges to cfsem's direct field at order 4.17 (1.4e-3 -> 7.6e-5 as dr halves);
  Delta* of cfsem's vacuum flux is zero to 3.9e-8 of Psi/dr^2.
- LU factorization at 161 x 161: 0.34 s (order 2), 1.0 s (order 4); reused across solves.

## Free-boundary equilibrium (`src/equilibrium.rs`)

Picard iteration on a fixed grid: Psi = Psi_coils (cfsem filaments, coils outside the grid) +
Psi_plasma. Each step finds the magnetic axis (the local flux maximum nearest the last one), the
boundary flux (the first limiter point touched, or an X-point inside the vessel if one sits
between), the plasma nodes (below that flux, inside the vessel, connected to the axis), and the
current J = lambda (beta0 R/R0 + (1 - beta0) R0/R)(1 - psiN)^2 scaled to Ip (the Jeon / FreeGS
profile with alpha1 = 1, alpha2 = 2, beta0 = 0.3). The plasma's own flux on the grid edge comes
from cfsem's filament flux; relaxation 0.5, stop at a 1e-8 relative change.

Two things it needs that a textbook sketch leaves out:

- Radial position feedback. With fixed coil currents only one vertical field holds the plasma;
  below it the plasma drifts out, scrapes the limiter, shrinks, needs even more field
  (ln(8R0/a)) and is lost (the first runs did exactly that). As in a real tokamak, the vertical
  field coils are set every iteration so the flux is equal at two midplane points (an isoflux
  constraint, exact since flux is linear in current).
- Up-down symmetry. J is averaged with its mirror image, since the vertical position is not
  controlled and elongated plasmas are vertically unstable. Asymmetric cases would need a
  vertical position control coil.

## Shafranov check (2026-11-02)

A far coil pair (R 3 m, Z +-2 m; decay index 0.088 at R0) under radial feedback centres a
circular plasma at R0 = 1 m, Ip 200 kA; grid 121 x 121 (1 cm). Its field at R0 against
Bv = -(mu0 Ip / 4 pi R0)[ln(8 R0/a) + beta_p + li/2 - 3/2] with beta_p and li from the solution:

| a/R0 | Bv(R0) | Shafranov | diff | beta_p | li | iterations |
|---|---|---|---|---|---|---|
| 0.3 | -55.53 mT | -57.08 mT | -2.71% | 0.313 | 1.515 | 73 |
| 0.2 | -64.13 mT | -65.02 mT | -1.38% | 0.308 | 1.509 | 72 |
| 0.1 | -78.42 mT | -78.59 mT | -0.22% | 0.299 | 1.496 | 69 |

The difference shrinks with a/R0, as the formula's O(a/R0) error should. At a/R0 = 0.2 the
plasma's grid flux matches cfsem's Green's sum for the same current to 3.4e-7 at 12,205 vacuum
nodes, and Ip integrates to 200.0 kA. Each equilibrium takes 1-2 s.

## Shape vs PF currents (2026-11-02)

R0 1 m, Ip 200 kA, elliptical vessel (0.42 x 0.68 m); grid 131 x 181. Coils: vertical field
pair (R 2.0, Z +-1.2, on radial feedback at R0 -+ 0.35 m), elongation pair (R 1.0, Z +-1.05),
triangularity pair (R 0.3, Z +-0.55), central solenoid (off).

| Case | VF (feedback) | kappa | delta | beta_p | li | touches |
|---|---|---|---|---|---|---|
| baseline | -136.6 kA | 1.060 | +0.081 | 0.315 | 1.530 | inboard wall |
| elongation 40 kA | -154.0 kA | 1.152 | +0.107 | 0.315 | 1.532 | inboard wall |
| elongation 80 kA | -171.2 kA | 1.292 | +0.133 | 0.318 | 1.550 | inboard wall |
| elongation 120 kA | -188.9 kA | 1.678 | +0.166 | 0.353 | 1.721 | top wall (0.92, 0.67) |
| elongation 80, triangularity +60 kA | -168.9 kA | 1.324 | +0.159 | 0.319 | 1.558 | inboard wall |
| elongation 80, triangularity -60 kA | -173.5 kA | 1.264 | +0.107 | 0.318 | 1.546 | inboard wall |

- Elongation follows the elongation pair; triangularity moves with the inner pair's sign (+0.159
  and +0.107 around +0.133 at 80 kA elongation; -60 kA lowers it but does not make it negative).
  The feedback adds vertical field as the shaping coils pull. R_geo stays at 0.983-0.986 m.
- At 120 kA the plasma reaches the top of the vessel just inside the separatrix (X-points sit
  just outside the vessel), hence the pointed top and bottom and the jump in kappa.
- 40-64 iterations a case, warm-started from the previous one.

## Plan

1. `check` (done).
2. Free-boundary equilibrium (done): radial feedback, symmetric shapes, shape follows the PF
   currents.
3. 3D: TF coils as cfsem linear filaments, field-line tracing in the total field, Poincare
   sections showing nested surfaces, q from tracing vs the equilibrium.
4. Output: plots, and CortenForge meshes of the coils and flux surfaces.
