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

## Probe results

None yet.

## Friction log

### `cortenforge` (facade) and docs overall

- **works** (2026-11-02): builds alongside cfsem 14 in one workspace, no
  dependency conflicts.

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

### `sim-rl`

- **bug** (2026-11-02): `Cem::best_artifact` returns parameters that were never evaluated (the
  elite mean after the update; cortenforge-play logged the same at cem.rs:200). Here it costs real
  design quality: on the spacing-only problem the shipped design scored 5.072 mm against 5.493 mm
  for the best design CEM had actually evaluated (-7.7%), and it shipped below its best evaluated
  design in all 20 CEM runs across three `optimize` runs (by 0.05-7.7%). Users have to log every
  evaluation in the reward closure to recover the true best.

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
- **API** (2026-11-02): `PrintValidation::summary()` gives counts only ("Not printable: 8215
  critical issue(s)"); grouping `issues` by type and reading `location` was needed to tell
  real overhangs from mesh artifacts.

### Platform

None yet.
