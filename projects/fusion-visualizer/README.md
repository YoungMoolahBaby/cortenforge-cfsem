# Fusion magnet visualizer (warm-up 2)

Trace field lines around a small tokamak or stellarator and show the plasma equilibrium shape.

- cfsem: circular filaments for the tokamak coils (TF, PF, central solenoid), linear filaments
  for stellarator coils, Grad-Shafranov (`cfsem::physics::gradshafranov`) for the equilibrium.
- CortenForge: coil and vessel geometry, mesh export of the coils and flux surfaces.
- Done when: field lines close on nested flux surfaces, and the equilibrium shape responds to the
  PF currents.
