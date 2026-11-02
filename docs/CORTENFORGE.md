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
  MuJoCo episode. From the source; the coil-designer optimizer step will try it.

### Platform

None yet.
