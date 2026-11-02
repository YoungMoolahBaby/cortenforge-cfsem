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
maglev-train started (2026-11-02), as an ultra-clean hybrid-EMS hospital delivery cart: `check`
validates the force and inductance model (19 of 19) and shows why hybrid (permanent magnets
carry the load, coils stabilize); `stray` sizes 12 x 4 mm discs and maps the 0.5 mT line;
`corner` levitates one corner in sim-core with zero-power control (80 mW with a 1 um gap sensor);
`cart` carries an off-centre payload 1.2 m to a station with modal zero-power control, driven by
a linear motor that pushes on the lift magnets (`motor`; 0.310 N/A, 0.10% ripple); `parts` gives
the printable cart and track parts, the landing-skid trade-off and a parts list. Main 2 done in
simulation, sideways guidance included (a second coil layer; a 0.3 N shove at the cart's end stays
0.64 mm off the 2 mm lips); open items in its README.
wireless-charging started (2026-11-02), to power the maglev cart from the track: `check` validates
cfsem's straight-filament inductance and field (9 of 9); cfsem's linear self-inductance counts the
wire's internal inductance twice (+5% on a single ring, corrected). Next: track loop and pickup.
magnet-tracking waits for hardware.
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
