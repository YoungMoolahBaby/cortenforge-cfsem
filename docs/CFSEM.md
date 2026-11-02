# cfsem 14.0.0: end-user findings

Friction found while using the published crate as an end user. CortenForge's findings live in
`CORTENFORGE.md`. Log it while building, each entry starting with its type (bug, docs, API, perf,
setup, works) and the date, with file:line in the crate source
(`~/.cargo/registry/src/*/cfsem-14.0.0/`).

## cfsem 14.0.0

- works, 2026-11-02: `physics::circular_filament::flux_density_circular_filament_scalar` gives the
  Helmholtz centre field to 3e-12 of the textbook formula (coil-designer smoke test).
