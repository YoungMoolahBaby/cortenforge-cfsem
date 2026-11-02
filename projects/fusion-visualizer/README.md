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

## Plan

1. `check` (done).
2. Free-boundary equilibrium: PF coils and a central solenoid as cfsem filaments, plasma current
   from p'(Psi) and FF'(Psi) profiles, Picard iteration with the boundary flux from all currents;
   show the shape moving with the PF currents.
3. 3D: TF coils as cfsem linear filaments, field-line tracing in the total field, Poincare
   sections showing nested surfaces, q from tracing vs the equilibrium.
4. Output: plots, and CortenForge meshes of the coils and flux surfaces.
