# Coil designer (warm-up 1)

Optimize Helmholtz or Maxwell coils for a uniform field (Maxwell coils: a uniform gradient).
The gentlest way to learn cfsem.

- cfsem: circular filaments, `flux_density_circular_filament_scalar` and friends.
- CortenForge: optimize radius, spacing and turns (sim-opt), export coil formers to print.

Now: a smoke test, `cargo run --release -p coil-designer`. It checks the Helmholtz centre field
against the textbook formula and prints the uniformity along and across the axis.
