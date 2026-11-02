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
wire's internal inductance twice (+5% on a single ring, corrected); `link` sizes the track loop and a 10-turn
pickup: 5 W from 5.6 A of track current at 85 kHz, 78% coil to coil, within +-2% over the trip
and the cart's motion; `effects`: the steel rail acts mostly as a conductor (6.61 A for 5 W, 0.6-1.2 W
of rail heating), lift coils pick up 0.3-1.1 V at 85 kHz (anti-alias filters needed), and the
field is under the public level 9 cm from the conductors; on the cart (maglev-train), steady running
needs about 2.0-2.3 W of the 4.3-4.5 W supplied, lift-off and a payload settle draw a buffer down
21.4 J (a 1.4 F supercapacitor bank, no battery), and the pickup former and the hangers' power
tray are printable. Main 4's cart side done in simulation. Pad and receiver started: `pad` models
the Qi A11 pad and receiver example 1 (inductances against the spec: +11% / +17% / +47% with
ferrite as ideal planes, an upper bound) and maps the sweet spot (k within 90% only 6 mm from
centre); `layout` widens it on a free-space efficiency score: one 80 mm coil covers 4x the A11's
area, three overlapping round coils 4.6x (2x Qi A6), with sim-opt's annealers matching or beating a
refined grid. `former` prints it: a base and lid at the A11's winding pitch (1484 mm^2, -3.6% for
slack between turns), both printable. Main 4 done in simulation.
magnet-tracking started simulation-first (2026-11-02): `check` validates the forward models (6 of 6):
cfsem's dipole to 1.2e-15, the magnet as cfsem loops against Derby and Olbert's closed form to
7.1e-6, and the dipole's error least at L/D 0.87 (theory 0.866); `solve` recovers a pose (9
unknowns) exactly, from a cold start in 7 ms, and reaches the Cramer-Rao bound with noise; an L/D
0.867 magnet keeps the dipole's bias under half the noise 15 mm out. Next: sensor layouts.
Hardware later.
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
