# CortenForge 0.9.0: end-user findings

Feedback gathered while using the published crates as an outside user would
(crates.io sources and docs only), here while pairing CortenForge with cfsem
for magnetics projects. cfsem's own findings live in `CFSEM.md`.

How to read it:
- **Dates** are when a finding was logged (local time). Newest at the bottom
  of each section.
- **Probe results** record what a probe or project measured, with the command
  that reproduces it. Every number comes from a real run.
- **The friction log** is grouped by crate. Each entry starts with its type
  and date: **bug** (wrong results or a crash), **docs** (missing or
  misleading docs), **API** (works, but awkward or missing a piece),
  **perf**, **setup**, or **works** (worth keeping as is). Paths are
  file:line under `~/.cargo/registry/src/*/<crate>-0.9.0/`. Bugs and doc gaps
  are the upstream fix list.
- **Suggestions** (at the end) are the nice-to-haves: features these projects wanted that no
  crate has, with what we built or used instead. Same format as cortenforge-play and
  aerial-animal-rescue.

## Probe results

- **sim-core's fluid forces at low Reynolds number** (2026-11-09; `cargo run --release -p
  microswimmers -- check`, section 6; glycerol, `viscosity="1.41"`, `density="0"`, 1 us steps).
  A free 1 mm sphere with `fluidshape="ellipsoid"` settles to Stokes' drag 6 pi mu a and torque
  8 pi mu a^3 (1.6e-14, 4.7e-15). With the default body model (inertia box) it is 22.5% low on
  drag and 53.5% low on torque, which is what that model's formula gives for a sphere (its
  equivalent box is 1.549 a a side). A free helix of 32 capsules (R 1.5 mm, pitch 6 mm, wire
  0.3 mm, 2 turns) twisted about its axis doesn't advance at all: 5e-18 of its pitch a turn
  (ellipsoid), 2e-21 (inertia box). A regularized-Stokeslet model checked to 1.25% against exact
  spheroids gives 0.2455. Its turning drag is 2.61x the model's (ellipsoid) and 3.81x (inertia
  box). Fluid density 1261 adds an advance under 2e-7, of either sign.

## Friction log

### `cortenforge` (facade) and docs overall

- **works** (2026-11-02): builds alongside cfsem 14 in one workspace, no
  dependency conflicts.
- **setup** (2026-11-06): the facade takes cf-design with `default-features = false`
  (cortenforge-0.9.0/Cargo.toml:76-80), which drops cf-design's default `parallel-meshing`, so
  `cortenforge::cf_design::Solid::mesh_adaptive_par` (solid/query.rs:175, documented as
  "default-on") does not exist through the facade. Nothing in the facade turns it back on. Workaround: also depend on
  `cortenforge-design = { version = "=0.9.0", features = ["parallel-meshing"] }`; Cargo's
  feature unification then enables it for the facade's copy (projects/stellarator/Cargo.toml).

### `sim-core`

- **works** (2026-11-02): a one-body maglev corner (maglev-train `corner`): an MJCF string with
  a slide joint and a box geom of given `mass`, `load_model`, `make_data`, `qfrc_applied` set each
  step from a cfsem force table, `data.step`. The open-loop magnet instability grows at 62.02 /
  62.04 1/s against the predicted sqrt(k/m) = 62.18 (within 0.25%, the same at 50 and 25 us
  steps), and the coil-shorted case matches the root of the electromechanical cubic equally well.
  The full mode (open loop, three 4 s closed-loop runs at 25 us, a six-run noise sweep) takes
  4 s.
- **docs** (2026-11-02): for a state-dependent external force the natural first move,
  `qfrc_applied`, is held through the step (a zero-order hold; fine here at 25 us). The
  callbacks that evaluate forces inside the step, `Model::cb_passive` and user actuator dynamics
  `cb_act_dyn` (types/model.rs:1015-1026), which could also carry a coil current as actuator
  state, are only found by reading `Model`'s fields; no guide or example mentions them.
- **works** (2026-11-02): those callbacks do the job well. The four-corner cart (maglev-train
  `cart`) evaluates its magnet forces inside `cb_passive` (a thread-safe `Arc<dyn Fn>` on the
  model, capturing the force table and shared coil state), maps each corner force to joint
  forces with `mj_jac_site` and the thrust and payload with `mj_jac` at arbitrary points, and
  holds 0.6 um agreement with an independently solved equilibrium. Four, later six, joints in one
  body (up to three slides and three hinges) and sites from MJCF behave as in MuJoCo; the
  six-joint version adds sideways and yaw control through the same callback. A 9 s trip at 25 us
  steps runs in about 4 s (plus the cfsem table building).
- **docs** (2026-11-02): `Data::xfrc_applied` is documented only as "Applied spatial forces in
  world frame" (types/data.rs:204); the source projects it as [torque, force] at the body's centre
  of mass (forward/acceleration.rs:128-133). MuJoCo's own layout is, as far as I recall, [force,
  torque], so code ported from MuJoCo would apply forces as torques silently. Not run here (the
  cart applies its point forces through `cb_passive` and `mj_jac` instead); one doc line would
  settle it.
- **works** (2026-11-09): `xfrc_applied` is [torque, force], now confirmed by a run
  (microswimmers `check` 6): a free sphere given a force in entries 3-5 moves without turning
  (0 rad/s) and settles to Stokes' drag to 1.6e-14. The docs gap above stands.
- **API** (2026-11-09): the viscous part of both fluid models is the same in every direction.
  Per geom it is beta 3 pi eq_d v, with eq_d from the mean semi-axis (forward/passive.rs:237-243);
  per body it is 3 pi beta d v (passive.rs:103-112). The direction-dependent terms (blunt,
  slender, Kutta, Magnus) all scale with the fluid's density, which is the fast-flow regime. So a
  slender piece is dragged equally along and across itself, and a helix built of pieces can't
  swim in syrup or as a microrobot (low Reynolds number), where swimming comes entirely from the
  drag across a wire being about twice the drag along it. Summing pieces as if each were alone
  also overstates drag (2.61x for a helix's turning; see Probe results). These are MuJoCo's
  documented models, so this is not a bug, but swimmers, flagella, cilia and particles settling in
  thick fluids are out of reach. See Suggestions.
- **works** (2026-11-09): a Stokes-flow body through `xfrc_applied` (microswimmers `steer`).
  - The setup: a free body with `<inertial ... fullinertia>` at its centre of mass, and each
    step the drag of the velocity the step will end with, v' = (M + dt R)^-1 (M v + dt F), with
    F the other known forces and R the body's 6 x 6 resistance matrix in its own frame.
  - Free of contact, sim-core's free-joint update lands on v' to 2e-4 at 1 ms steps (4.7e-6 at
    20 us). `Data::qM` matches our mass matrix to 1.2e-20.
  - A 70 mg swimmer that settles in 0.1-0.55 ms runs at 1 ms steps (about 36x faster than real
    time with 143 contact spheres), with its speed within 0.06% of the 20 us result.
  - Linear velocity is in the world frame and angular velocity in the body frame, as in MuJoCo.
  - The limit: the drag ignores contact, so in contact friction is too weak at long steps. A
    swimmer pressed against a wall slides at 9.6 mm/s with 1 ms steps, 7.2 at 0.25 ms and 6.4
    with plain drag at 20 us. See Suggestions (Stokes-flow bodies).
- **bug** (2026-11-09): capsule-box contact pushes the wrong way. A free 1 kg capsule (5 cm
  radius) dropped flat on a box sinks to the box's centre and stays there, reporting depth =
  radius, whichever geom is listed first (`cargo run --release -p microswimmers -- steer
  contact`; 2 ms steps, default contact).
  - The likely cause: `collide_capsule_box` sets `normal_sign = +1` when the capsule is geom1
    (collision/pair_cylinder.rs:305). But its normal points from the box to the capsule, and
    sphere-box negates exactly that case so the normal runs geom1 to geom2
    (pair_convex.rs:360-364).
  - In the swimmer course, capsules were pulled 2-4 mm into the glass's boxes, and once shot
    through a fin at 1 m/s.
  - Sphere-box contact is fine.
- **bug** (2026-11-09): some cylinder contacts are never found, so bodies fall straight through
  (the same probe).
  - A cylinder dropped on a box, on a cylinder's flat end or on a cylinder's side: no contacts,
    and it ends 30-44 m below.
  - A capsule on a cylinder's flat end does the same.
  - Cylinder pairs go to the general GJK/EPA path (collision/narrow.rs:170), and capsule-cylinder
    falls back to it when its own routine finds nothing (narrow.rs:156-168).
  - A capsule on a cylinder's side, and anything on a plane, are fine.
- **works** (2026-11-09): sphere contacts with planes, boxes and cylinders, at an everyday scale
  and at the swimmer's (70 mg, 0.25-1.45 mm spheres, default `solref` 0.02): pressed into a wall
  by its own thrust, the swimmer goes 0.006-0.07 mm in. The workaround for the two bugs above is
  a chain of spheres along each capsule.

### `sim-opt`

- **API** (2026-11-02): the optimizers (`Sa`, `RicherSa`, `Pt`) only run through
  sim-ml-chassis's `Algorithm::train` on a `VecEnv`, and a `VecEnv` is built from a sim-core
  `Model` (sim-ml-chassis vec_env.rs:91, 118). There is no black-box `minimize(f, x0)` for
  design parameters, so a static problem like coil spacing has to be dressed up as a one-step
  MuJoCo episode. Tried in coil-designer `optimize`; see below.
- **works** (2026-11-02): the dressing-up works and costs nothing at run time. A dummy MJCF
  model (one body, two motors with `gear="0"`), `ObservationSpace` on time with the policy's obs
  scale 0, `ActionSpace` on ctrl 0..2, a reward closure that decodes `ctrl` and scores the coil,
  `done` always true: about 50 lines (coil-designer src/optimize.rs), compiled first try. 400
  evaluations take 3.2-3.5 s, the same as 400 bare scores (8.4 ms each), so the env adds no
  measurable overhead.
- **API** (2026-11-02): what the dressing-up costs in the design itself. Parameters arrive
  through `LinearPolicy`'s tanh and an f32 action boundary (sim-ml-chassis rollout.rs:115-121,
  space.rs:695), so a bound is only reached as |b| -> infinity and the resolution
  is about 1e-7 of the range. The policy needs an observation, so two dead weights ride along in
  every proposal. There are no integer or bounded parameters; the reward has to decode, round
  and penalize infeasible designs itself. A plain parameter-vector policy (identity, clamped to
  ctrlrange) would remove most of this.
- **works** (2026-11-02): on one continuous parameter (pair spacing at a fixed 10 x 10 pack, a
  sharp ridge, 200 evaluations, 3 seeds) `RicherSa` finds the peak every time, 0.10-0.14% *above*
  a 5e-5 R brute-force grid; `Sa` median -0.23% (worst -2.6%), `Pt` median -0.50% (worst -1.5%).
  The 1/5 rule is what makes the difference on a ridge this narrow.
- **API** (2026-11-02): on two coupled parameters (spacing and an integer layer count, 1-25)
  all four optimizers (`Sa`, `RicherSa`, `Pt`, sim-rl `Cem`) end 24-34% below the brute-force
  optimum, and 4x the budget (400 -> 1600 evaluations) barely moves them: no seed found the
  winning single-layer design. Each layer count has its own narrow best spacing (0.97-1.57 R),
  so a step in layers without the matching step in spacing falls off the ridge and is rejected.
  That is the problem's shape, not a bug, but the docs pitch these as general gradient-free
  optimizers. A grid over the integer with a 1-D search inside (10,800 evaluations, 90 s here)
  wins easily; a note that they are local, isotropic and continuous-only would set expectations.
- **works** (2026-11-02): on a smooth score in 3-4 parameters the annealers are as good as an
  exhaustive search when a run doesn't stick (next entry). Wireless-charging `layout` sizes a
  charging pad (inner and outer radius, turn count, array pitch; the score is a logistic-weighted
  area, so it has no cliffs): with 1600 evaluations (1.7-2.2 s) `RicherSa` reaches a 0.1 mm
  refined grid's optimum for one coil (1301.9, all three seeds within 0.04%), and `RicherSa` and
  `Sa` pass it for a three-coil array (1522.5 / 1523.0 vs 1518.6), finding a second optimum
  outside the refinement's window. Both also beat the plain 1 mm grid I first used as the
  reference (1316 vs 1276 mm^2). The contrast with coil-designer is the landscape (smooth here;
  one narrow ridge per integer there), which the docs could say.
- **API** (2026-11-02): `Sa` and `RicherSa` have no restart, and either can stick: on the
  three-coil array one `RicherSa` seed in three settled at 35.6 (98% below the others) at both
  400 and 1600 evaluations, and one `Sa` seed at 48.3 at 400; in an earlier version of the same
  problem (a looser constraint) it was `Sa`'s seed that stuck, at 35.8. Each sits in a region the
  design decoder scores near zero. A restart on stagnation, or multiple starts, would save users
  running seeds by hand.
- **works / API** (2026-11-09): `RicherSa` on a 3-parameter swimmer design (microswimmers
  `shape`: helix radius, pitch angle and wire, scored by the climb at 2 mT; 60 evaluations, seed
  1, one run for each of 2, 3 and 4 turns; a 4 x 4 x 3 grid as the reference). It beat the grid
  at 2 turns (+2.10 against +1.86 mm/s) and 3 turns (+0.76 against +0.66). At 4 turns it stuck
  at a wide, steep helix (-0.74, where the grid found +0.96), the same kind of stuck seed as the
  entry above. The dummy-model setup carried over unchanged from coil-designer. A cheap grid
  alongside, keeping the better result, was the guard.

### `sim-rl`

- **bug** (2026-11-02): `Cem::best_artifact` returns parameters that were never evaluated (the
  elite mean after the update; cortenforge-play logged the same at cem.rs:200). Here it costs real
  design quality: on the spacing-only problem the shipped design scored 5.072 mm against 5.493 mm
  for the best design CEM had actually evaluated (-7.7%), and it shipped below its best evaluated
  design in all 20 CEM runs across three `optimize` runs (by 0.05-7.7%). Users have to log every
  evaluation in the reward closure to recover the true best. Seen again in wireless-charging
  `layout` (2026-11-02): all 12 CEM runs shipped below their best evaluated design (by
  0.1-1.2%), and at 1600 evaluations ended 0.3-3.1% under a refined grid where `RicherSa`
  matched it.
- **works** (2026-11-06): CEM on a 45-parameter design (stellarator `optimize`: three free coils'
  Fourier displacements, scored by field-line tracing at about 10 s an evaluation). The
  dummy-MJCF wrapper from coil-designer scaled from 2 to 45 zero-gear motors with no change but
  generating the MJCF. 24 generations x 24 candidates (noise 0.2, decay 0.93): the population mean
  rose from -6.0 to 72.1, the best evaluated design (72.50, the score's ceiling) came at
  evaluation 262 of 577, and 28 candidates (5%) were degenerate (no magnetic axis). Every
  evaluation logged, so `best_artifact` (entry above) was not needed.

### `cf-design`

- **works** (2026-11-02): code-first CSG is pleasant for a coil former: cylinders, a flipped
  cone, `union`/`subtract`/`intersect`/`translate`/`rotate`, about 80 lines for two printable
  parts (coil-designer src/design.rs). The exact SDF doubles as a measuring tool: probing the
  solid at the floor, winding surface, bore and mid-plane faces reads 0.0000 mm, and evaluating
  it at every mesh vertex measures mesh error directly.
- **bug** (2026-11-02): the sharp-feature meshers misplace vertices at concave corners.
  `mesh_dc` and `mesh_adaptive` (solid/query.rs:113, 141) on the one-piece former (cylinders
  unioned and subtracted, 0.5 mm cells) leave 1% of vertices more than a quarter cell off the
  surface, worst 0.48 / 0.46 mm, all where a flange meets the tube, pulled into the solid. The
  folds read as 2,573 / 3,716 "critical" thin walls of 0.03 mm in mesh-printability. Marching
  cubes (`mesh`, :82) on the same part: worst 0.124 mm (a quarter cell, corner rounding), no
  thin walls. The docs pitch DC as "preserving sharp features"; on concave corners it does
  worse than MC.
- **bug** (2026-11-02): `mesh` ("watertight, manifold", solid/query.rs:75) is both by edge
  count, but emits coincident vertices under different indices (914 pairs on a former half at
  0.5 mm cells, 1,489 at 0.25), giving 1,828 / 2,978 zero-area faces, plus sub-micron slivers
  where the surface grazes grid corners (an edge of 0.59 um). mesh-printability reads these as
  100+ self-intersecting pairs. At 0.5 mm cells `mesh_repair::weld_vertices` at 1 um clears
  them with no face collapsing. At 0.25 mm it leaves 13 zero-area faces and 34 pairs of slivers
  2 um apart, all on the 40 deg support cone where it crosses grid planes (z = 10.75, 11.25,
  ...), and the remaining "critical" overhangs (10) and thin walls (6) are 0.0 mm^2 single faces
  at the same planes. No weld setting tried makes `is_printable()` true for the half, though
  nothing flagged is a real print problem. A weld and sliver cleanup inside `mesh` would save
  every user this detour.
- **perf** (2026-11-02): `mesh` (single-threaded) on the former half, 114 mm across: 11 s at
  0.5 mm cells, 85 s at 0.25 mm (1.06M faces, 53 MB STL). `mesh_adaptive` at 0.25 mm took 24.5 s
  for 2.4M faces on the one-piece former.
- **works** (2026-11-02): coil solids for a tokamak: rectangular-section rings (cylinder minus
  cylinder) and a torus rotated into a poloidal plane, meshed by marching cubes and welded: all
  watertight and manifold, volumes within -1.7e-4 to -2.3e-3 of the exact formulas at 5 mm cells
  (-7.6e-4 for the torus at 2.5 mm). The 4.15 m vertical-field ring took 3.1 s (671k faces). One
  torus meshed once and rotated 16 times gave the TF set.
- **API** (2026-11-02): no tube or annulus primitive, and `Solid::cylinder` panics on radius 0
  (solid/primitives.rs:58), so a helper built on subtract needs a special case at r = 0.
  `Solid::cone` opens downward from its apex (:153); an upward support cone needs a rotation.
- **API, recommendation** (2026-11-02): add a native revolve, e.g. `Solid::revolve(profile)` for
  a closed (r, z) profile (a 2D SDF or a polyline) swept about Z. In an SDF kernel it is one of
  the simplest operations (evaluate the profile at (sqrt(x^2 + y^2), z); exact wherever the
  profile's SDF is) and it composes with the CSG, shell, meshing and printability tools already
  here. Most parts in this repo were solids of revolution and each needed a workaround: the coil
  former (cylinder-minus-cylinder rings and a flipped, clipped cone for a 40 deg support), the
  tokamak's PF rings (the same), the vessel and flux surfaces (hand-built meshes, since their
  profiles are curves). Bobbins, pulleys, nozzles, flanges and spools are the same class.
  What exists today:
  - `Solid::loft(stations)` (solid/primitives.rs:267-281) is a partial revolve: a radius r(z)
    swept about Z. It is always solid to the axis (no bore, so no ring, torus or vessel), its
    radius is Catmull-Rom smoothed between stations (steps round off), and its SDF is approximate
    where the radius changes. Its docs don't present it as a revolve.
  - `Solid::user_fn` (solid/ops.rs:499-516) can carry a hand-written revolve, and it works: the
    vessel profile (160-point polygon) revolved this way meshes watertight and manifold at 10 mm
    cells, its volume -6.7e-5 from Pappus (fusion-visualizer `mesh`). But it took 16.4 s, against
    under a second for the same vessel as a hand-built revolved mesh: a `user_fn` has no interval
    function, so the mesher cannot prune empty space, and each evaluation scans every profile
    segment. A native revolve could supply both.
- **bug** (2026-11-02): `mesh` output depends on how the faces fall on its grid. The maglev
  hangers with the power tray (boxes unioned and subtracted, printed on end), meshed at 0.4 mm
  and welded at 1 um: 12 and 18 self-intersecting triangle pairs, a 0.27 mm sliver flagged as a
  thin wall (0.0 mm^2) and, for the station hanger, 4060 non-manifold edges. The same solids at
  0.37 mm: watertight, no self-intersections, no thin walls (maglev-train `parts`). The grid
  starts one cell outside the bounds (solid/query.rs:93-96), so faces a whole number of cells
  from the bounds land on grid nodes; that is the suspected cause, not confirmed. Users can't
  see this without a validator, and the fix (another cell size) is found by trial.
- **works** (2026-11-02): a pickup former for an 85 kHz coil (a plate minus a rectangular ring
  channel and two lead holes): probes at the channel floor, both channel edges and a hole wall
  read 0.000 mm; printable, 1.18M faces in 12.4 s at 0.2 mm.
- **API** (2026-11-06): `Solid::pipe_spline` (solid/primitives.rs:243) is open-curve only
  (Catmull-Rom with clamped end tangents, cortenforge-geometry curve.rs:224-233); there is no
  closed option, though a closed loop is the commonest swept shape for coils. Workaround (stellarator
  `mesh`): wrap the control points as [p_(n-1), p_0, ..., p_(n-1), p_0, p_1], so every true span
  has its real neighbours (the two extra end spans repeat true spans). The SDF then matches the
  true tube distance to 0.060 mm at 64 points and 0.0059 mm at 256 (4 mm tube on a 460 mm
  ellipse, 20,000 random points within 6 mm; `cargo run --release -p stellarator -- sdf-probe`).
- **perf** (2026-11-06): `mesh` (uniform marching cubes) is unusable for a long thin sweep: one
  4 mm-radius coil 460 mm long took 305 s at 1 mm cells (32 control points). Its per-cell
  "interval pruning" (mesher.rs:83-87) calls `evaluate_interval`, which for `PipeSpline` is
  `lipschitz_bound` (interval.rs:62): a full SDF evaluation at the cell centre. So every cell of
  the coil's bounding box costs at least one full spline distance (every span, 17 samples plus
  Newton). The octree `mesh_adaptive_par` does the same coil in 1.7 s; at 64 points, 3.9 s against
  20.2 s for the sequential `mesh_adaptive` (identical output). A note on `mesh` pointing sweeps to the octree mesher would save
  the trial.
- **bug** (2026-11-06): the octree dual-contouring meshers put vertices well off a smooth surface.
  On the same tube (no sharp features), `mesh_adaptive_par` and `mesh_adaptive` give identical
  meshes, watertight and manifold, median vertex 0.0075 mm off the tube at 0.5 mm cells but 9% of
  vertices more than 0.1 mm off and the worst 0.444 mm (0.9 of a cell); at 1 mm cells, 0.034
  median and 1.025 mm worst. The solid's own `evaluate` reads up to 0.444 mm at those vertices,
  so the SDF is right and the vertex placement is wrong (`stellarator -- mesh-probe`, `sdf-probe`).
  Same family as the concave-corner entry above, here on a smooth surface. Workaround: two Newton
  steps per vertex along the SDF gradient (`evaluate` and central differences) put every vertex
  within 0.015 mm of the true tube and keep the mesh watertight and manifold; 1 mm cells then
  give 31-35k faces a coil. The remaining volume deficit is the facets cutting chords across the
  round tube: -0.87% at 1 mm, -0.20% at 0.5 mm (second order).
- **docs** (2026-11-06): the octree meshes come out clean (0 boundary, 0 non-manifold, 0
  duplicate faces), but the `weld_vertices` at 1 um that marching-cubes output needs (entry
  above) makes them non-manifold (914 edges at 1 mm, 1,506 at 0.5 mm, with 431 / 763 duplicate
  faces). Neither mesher's docs say which post-processing it expects.
- **bug** (2026-11-06): `mesh_adaptive_par`'s topology depends on where a part sits. The
  optimized stellarator coils of one field period mesh watertight and manifold, but the same
  coils turned by 120 or 240 degrees (exact copies) gave 8-10 non-manifold edges (and 2
  degenerate faces) on 4 of 12, at 1 mm cells, 128 control points (stellarator `best`, before the
  fix). Snapping moves vertices only, so the edges come from the mesher. Like the uniform
  mesher's grid-placement entry above. Workaround: mesh one period and rotate the meshes (also 3x
  faster: 57 s for 18 coils).
- **bug** (2026-11-09): `mesh_adaptive_par` folds triangles on a thin, curved `pipe_spline`
  (`cargo run --release -p microswimmers -- shape probe`). The case is a 2-turn helix of wire
  radius 0.254 mm (helix radius 1.56 mm, pitch 12.2 mm) at 0.051 mm cells (5 per wire radius).
  The mesh is watertight and manifold, but `validate_for_printing` (SLA) reports 101 overlapping
  triangle pairs and 4,681 thin-wall regions of 0.000 mm. A straight `pipe_spline` or a capsule
  of the same radius at the same cells has none (thinnest 0.117 mm raw, 0.279 mm snapped).
  Snapping the vertices onto the surface doesn't help (102 pairs, 4,915 regions). The uniform
  `mesh` on the same helix, welded at 10 nm, has no thin walls and 26 overlap flags, but takes
  695 s against 6 s. On the whole swimmer (helix, strut, head, pocket) the octree mesher gives
  102 pairs and 4,815 thin regions in 155 s; uniform marching cubes give 31 flags and 0 thin
  regions (thinnest 0.397 mm, the 0.4 mm wall as designed) in 918 s. The flags left on the
  welded marching-cubes meshes are probably the touching slivers of the `mesh` self-intersection
  entry below; I didn't inspect them. Workaround: the uniform mesher for the printed part
  (microswimmers `shape part`).

### `mesh` (printability, repair, io)

- **works** (2026-11-02): a tokamak as meshes (fusion-visualizer `mesh`): IndexedMesh built
  directly for revolved flux surfaces, cf-design solids for the coils, `validate_mesh` and
  `weld_vertices` from mesh-repair, `save_mesh` picking STL from the extension. All 14 parts come
  out watertight and manifold. `signed_volume` and `surface_area` make good checks: revolved
  surfaces match Pappus to exactly the faceting term (-1.00e-4 at 256 sides, (2 pi / n)^2 / 6),
  and the plasma surface mesh matches the equilibrium's own volume sum to 6.2e-4.
- **API** (2026-11-02): no mesh-level lathe (revolve a polyline into a mesh). `mesh-loft` sounds
  like one, but it joins two contact patches into a bushing (mesh-loft lib.rs:1-30). For profiles
  that already exist as polylines (contours, measured outlines) it is the natural tool, with exact
  facets; here it was 40 lines over `IndexedMesh`. The solid-modelling side is under `cf-design`
  above (recommendation: a native revolve).
- **docs** (2026-11-02): volume and area live on cf-geometry's `IndexedMesh`
  (`signed_volume`, `surface_area`; cortenforge-geometry mesh.rs:202, 255), not in `mesh-measure`,
  the crate named for it (dimensions, OBB, cross-sections, distances only). Found by grepping.

- **works** (2026-11-02): the overhang check caught the real flaw in the first former: printed
  upright, three 12 mm annular channel ceilings (3,820 mm^2 each) at 90 deg. That drove the
  redesign (two halves printed on their mid-plane faces, support cones under the floors, flat
  washers), after which the washers are "ready for printing" and the half has only sub-mm^2
  artifacts left (below).
- **bug** (2026-11-02): self-intersection false positives. The check skips pairs that share an
  edge by index and says "vertex-only contact is not an intersection ... and is not flagged"
  (mesh-printability validation/mod.rs:2329-2331), but all 100 reported pairs on an unwelded
  marching-cubes mesh touch at coincident positions under different indices, and after a 10 nm
  weld the rest are coplanar sub-micron slivers touching at a point. Each run also reports a
  different first pair when truncated at 100. Either weld before the check or test touching
  geometrically.
- **docs** (2026-11-02): a flat annular overhang is reported as a `LongBridge` "spanning 113.7
  mm" (the ring's diameter; the span is "max of the two perpendicular-plane bbox extents",
  regions.rs:101). It is a 12 mm cantilever, not a bridge; the message points at the wrong fix.
- **API** (2026-11-02): `find_optimal_orientation` (orientation.rs:85) tilted the one-piece
  former 45 deg (all faces then at or under 45 deg, but resting on a flange rim, no bed
  contact), as cortenforge-play found for its parts. On the redesigned parts it keeps them as
  modelled.
- **works** (2026-11-02): the maglev cart's deck and hull pass the FDM check as modelled, and a
  constant cross-section printed on end (the rail hanger) is a reliable way to avoid overhangs.
  The only remaining flags are marching-cubes-blunted apexes of diamond holes (about 0.4 mm wide,
  3.6 and 4.3 mm^2 regions at the hole tops), the same pattern as the coil former's dowel tips; a
  thin wall left by a mis-placed hole (0.25 mm) showed up as a non-manifold mesh, which the check
  reported (maglev-train `parts`).
- **API** (2026-11-02): `PrintValidation::summary()` gives counts only ("Not printable: 8215
  critical issue(s)"); grouping `issues` by type and reading `location` was needed to tell
  real overhangs from mesh artifacts.
- **works** (2026-11-02): `validate_for_printing` caught a real flaw, not a false positive this
  time. Wireless-charging `former` cuts two coil channels (annuli) into a plate that nearly touch
  at the plate's centre; the check flagged two "ThinWall" regions of about 400 faces (0.003 and
  0.005 mm at the cusp, at x = 0, y = +-3.4 mm, located correctly). The wall between the channels
  is under 1 mm for |y| < 5.2 mm; joining the channels there makes the part printable. The
  location field is what made it quick to act on.
- **works** (2026-11-09): with `PrinterConfig::sla_default` on the microswimmer (microswimmers
  `shape part`, marching cubes welded at 10 nm), every thin-wall flag landed on the one wall
  designed at the preset's 0.4 mm minimum (0.397-0.398 mm, faceting), as warnings, and nothing
  else. The preset's 0.4 mm is what set the head's size (2.9 mm across a 2 mm magnet).
- **bug** (2026-11-09): the self-intersection count isn't stable. The same part with design
  numbers differing only in the fifth digit gave 31 and 92 overlapping pairs (welded marching
  cubes, 0.051 mm cells), which fits the touching-sliver false positives above rather than real
  crossings. The issue carries no location (`location: None`), so the pairs can't be found from
  the report.

### Platform

None yet.

## Suggestions

- **Linear-elastic structures** (2026-11-06): small-strain linear elasticity for metal parts
  (isotropic E and nu), a static solve and natural modes, body loads (gravity, the centrifugal
  rho w^2 r of a spinning structure, magnetic body forces from cfsem's `body_force_density`
  functions), stress out (von Mises, safety factor), on cf-design solids or their meshes.
  sim-soft is the only FEM here and is hyperelastic soft-body dynamics: `NeoHookean` and `Yeoh`
  (sim-soft lib.rs:1-14), silicone tables, backward Euler. Neo-Hookean reduces to linear
  elasticity at small strain, so sim-soft might take steel at 200 GPa, untested (it is tuned for
  materials 10^5 times softer). Wanted for: stellarator coil supports under magnetic loads, and a
  spinning space habitat (rings, spokes, hub). Meanwhile: hand formulas (thin rotating rings and
  pressurized shells), then the `pmsim` crate (solids and structures, pure-Rust default build) as
  an outside check, not tried yet. The mothership's active shield (cortenforge-cfsem active-shield)
  needs it too: its carbon-fibre structure holds 3-12e14 J of magnetic pressure, and its mass
  scales with the virial factor Q (3 with ties, ~10 as hoops), which only a stress solution of
  the real windings and ties can pin down.

- **A native revolve** (2026-11-02): `Solid::revolve(profile)`; see the `cf-design` entry above.
  Built instead: a revolve on `Solid::user_fn` (fusion-visualizer `mesh`, 16 s to mesh) and
  hand-built revolved meshes.
- **Closed sweeps** (2026-11-06): `pipe_spline` (and a sweep with a non-round section) along a
  closed curve; coils are closed loops. Built instead: wrapped control points (stellarator `mesh`).
- **A plain parameter optimizer** (2026-11-02): `minimize(f, x0, bounds)` for design
  parameters, with bounded and integer parameters, restarts on stagnation, a population
  evaluated in parallel, and CEM returning its best evaluated design. Built instead: a dummy MJCF
  model with one zero-gear motor per parameter (coil-designer `optimize`, stellarator
  `optimize` with 45 parameters) and a log of every evaluation to recover the true best.
- **Mesh post-processing for the octree mesher** (2026-11-06): snap vertices onto the SDF's zero
  set, as an option of `mesh_adaptive_par`. Built instead: two Newton steps a vertex along the
  SDF gradient (stellarator `mesh`), 0.44 mm worst error down to 0.015 mm.
- **Stokes-flow bodies** (2026-11-09): low Reynolds number (microrobots, swimmers in syrup,
  particles settling). Either (a) slender-body viscous drag per geom, with different
  coefficients along and across a capsule set from its aspect ratio as in resistive force
  theory, or better (b) a body whose fluid force is a given 6 x 6 resistance matrix times its
  velocity relative to the fluid, with the matching velocity derivative so implicit integrators
  stay stable when a light body sits in a thick fluid. (b) also captures shielding between a
  body's own parts, which (a) can't (see the `sim-core` API entry).
  - Built instead: a regularized-Stokeslet model in plain Rust (microswimmers `check`, within
    1.25% of exact spheroids).
  - Its resistance matrix goes in through `xfrc_applied` as a linearly implicit drag
    (microswimmers `steer`; see the `sim-core` works entry). That is exact and stable free of
    contact, but the contact solve doesn't see the drag, so friction is too weak at long steps.
  - (b) done natively would add dt R to the mass matrix, as the implicit integrators already
    build M + h D for joint damping (forward/acceleration.rs:149, `M_impl`). Used in the
    constraint solve as well, that would make contact consistent too.
