# Findings

Friction found while using the published crates as an end user. Log it while building, under
its crate, each entry starting with its type (bug, docs, API, perf, setup, works) and the date,
with file:line in the crate source (`~/.cargo/registry/src/*/<crate>-<version>/`).

## cfsem 14.0.0

- works, 2026-11-02: `physics::circular_filament::flux_density_circular_filament_scalar` gives the
  Helmholtz centre field to 3e-12 of the textbook formula (coil-designer smoke test).

## cortenforge 0.9.0

- works, 2026-11-02: builds alongside cfsem 14 in one workspace, no dependency conflicts.
