# Plan

Status: coil-designer started (2026-11-02): `check` validates cfsem and the turn-by-turn coil model
(58 of 58 pass); `pack` sizes a real Helmholtz pair; `optimize` compares sim-opt/sim-rl with a
grid; `design` picks the least-power pair meeting a robust spec and builds its printable
former (with lead slots). Next: print, wind and measure (see its README).
fusion-visualizer started (2026-11-02): `check` validates cfsem's Grad-Shafranov operator and flux
conventions (16 of 16); a free-boundary solver matches Shafranov's vertical field (-0.22% at
a/R0 0.1) and its plasma shape follows the PF currents; field lines traced through 16 TF coils
lie on nested flux surfaces with q matching the equilibrium; coils, vessel and flux surfaces
export as checked CortenForge meshes. Warm-up 2 done; next in the plan: magnet-tracking.
The other projects are a README each until they start.

## Order

1. coil-designer: learn cfsem's filament API and check it against textbook fields.
2. fusion-visualizer: field lines and Grad-Shafranov.
3. magnet-tracking: the first project with hardware; its solver feeds colon-probe.
4. maglev-train or wireless-charging, then soft-robots.
5. colon-probe: uses tracking (localization) and charging (capsule power).
6. Extras as time allows.

## How the projects connect

- magnet-tracking -> colon-probe (localization)
- wireless-charging -> colon-probe phase 2 (capsule power)
- coil-designer -> every project that drives coils (maglev, conveyor, microswimmers, colon-probe)
- magnet-arrays -> maglev-train, planar-conveyor (magnet models)

## Starting a project

Add `projects/<name>/Cargo.toml` and `src/`, add the path to `members` in the root `Cargo.toml`,
and take cortenforge/cfsem with `.workspace = true`.
