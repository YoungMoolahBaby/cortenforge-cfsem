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

### Platform

None yet.
