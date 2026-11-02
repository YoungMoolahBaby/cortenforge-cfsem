# cortenforge-cfsem

Example projects joining CortenForge with cfsem. Read first:
- `README.md`: the projects.
- `docs/PLAN.md`: order, status, how the projects connect. **Start here.**
- `docs/CORTENFORGE.md`: end-user findings for CortenForge (improvement notes). Add to it while building.
- `docs/CFSEM.md`: end-user findings for cfsem. Add to it while building.

## Rules

- Use CortenForge and cfsem only as an end user would: the published crates and their
  sources/docs under `~/.cargo/registry/src/*/`. Don't use their git repos.
- CortenForge and cfsem are the focus; other crates only render, display and wire things up.
- Check physics against known results (textbook formulas, cfsem's own tests) before building on it.
  Every number in docs comes from a real run.
- Log friction in `docs/CORTENFORGE.md` (CortenForge, under its sub-crate) or `docs/CFSEM.md`
  (cfsem): type (bug, docs, API, perf, setup, works), date, file:line.

## Build notes

- Windows 11, RTX 4070 Ti, Ryzen 5 5600. Always run in release: `cargo run --release -p <project>`.
  The first release build takes ~3-4 minutes.
