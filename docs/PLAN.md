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
magnet-tracking started simulation-first (2026-11-02): `check` validates the forward models (5 of 5):
cfsem's dipole to 1.2e-15, the magnet as cfsem loops against Derby and Olbert's closed form to
7.1e-6, and the dipole's error least at L/D 0.87 (theory 0.866); `solve` recovers a pose (9
unknowns) exactly, from a cold start in 7 ms, and reaches the Cramer-Rao bound with noise; an L/D
0.867 magnet keeps the dipole's bias under half the noise 15 mm out; `layout` picks a 4 x 4
MMC5983MA board (under 1 mm to 120 mm) and shows calibration (gains to ~0.1%, positions under
0.1 mm), not sensor noise, limits a real board; `calib` recovers a board's 240 unknowns from a magnet
in a jig, and a joint solve (bundle adjustment) makes even a printed jig (0.2 mm, 0.5 deg) good to
32-70 um; `path` tracks a moving magnet at the bound frame by frame, and a solver that models the
motion removes the bias of reading sensors in turn (to 10 um at 200 mm/s). Simulation side done;
next with hardware.
Hardware later.
stellarator started (2026-11-06), a tabletop vacuum-field stellarator: `check` validates the
Fourier coils, cfsem's straight-filament field for them and an RK4 field-line tracer (21 of 21);
`scan` shows planar rotating-ellipse coils give iota 0.06-0.10 and picks 18 coils (3 periods,
100 x 40 mm ellipses, a 20 mm helical shift); `trace` finds the axis as the period map's elliptic
fixed point and nested surfaces out to about 52 mm, iota 0.105 on axis to 0.082 at the edge
(ragged at iota = 1/12); `mesh` sweeps the coils with cf-design's `pipe_spline`, meshes them with
the octree mesher and snaps the vertices onto the surface (within 0.015 mm, watertight and
manifold, STL). Coil search run 1 (`optimize`, sim-rl CEM on 45 Fourier coefficients of three
free coils, the rest by symmetry; 577 evaluations, 94 min): iota 0.105 -> 0.269 on the axis
(rising to 0.293 at the edge), good region 52.5 -> 65 mm, coils bent up to 14.5 mm out of plane;
the score saturated, so next is a score that can't, then coil forces, formers and a build.
active-shield started (2026-11-06), a magnetic radiation shield for the mothership
(C:\dev\mothership): a confined toroidal field between two windings round its 70 m tube. `check`
validates cfsem's oriented-loop field for 2,000 loops, stored energy, magnetic pressure and
centering forces, a Boris tracer and a free-space cosmic-ray dose model (15 of 15);
`transmission` backtraces cosmic rays (half get through at ~1.1 x c B d; thicker, weaker shells
are sharper and lighter); `size` finds a magnet + water hybrid ~3x lighter than water alone (5.6-
5.8 Mt against 18.6 Mt for 10 mSv/yr), at the cost of ~131 MW of cryocoolers and 3-12e14 J
stored; `failsafe` sizes the mothership's design: 4.72 t/m^2 of water (50 mSv/yr if the magnet
fails) plus a 50 T m magnet run warm at ~65 K, 11.4 Mt, 55 MW, ~18 kt of TNT stored. Next:
transport for the secondaries, the support structure, quench protection.
The other projects are a README each until they start.

## Order

1. coil-designer: learn cfsem's filament API and check it against textbook fields.
2. fusion-visualizer: field lines and Grad-Shafranov.
3. magnet-tracking: the first project with hardware; its solver feeds colon-probe.
4. maglev-train or wireless-charging, then soft-robots.
5. colon-probe: uses tracking (localization) and charging (capsule power).
6. Extras as time allows.
7. stellarator: fusion-visualizer's field lines in 3D coils; its coil search uses what
   coil-designer learned about sim-opt and sim-rl.
8. active-shield: the mothership's magnetic shield (oriented loops, a particle tracer, dose);
   feeds its `habitat.shield.active_*` parameters.

## How the projects connect

- magnet-tracking -> colon-probe (localization)
- wireless-charging -> colon-probe phase 2 (capsule power)
- coil-designer -> every project that drives coils (maglev, conveyor, microswimmers, colon-probe)
- magnet-arrays -> maglev-train, planar-conveyor (magnet models)

## Starting a project

Add `projects/<name>/Cargo.toml` and `src/`, add the path to `members` in the root `Cargo.toml`,
and take cortenforge/cfsem with `.workspace = true`.
