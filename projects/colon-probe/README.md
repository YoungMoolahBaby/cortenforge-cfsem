# Magnetically steered colonoscopy probe (main 3)

Steer a probe tip through a colon model with an external magnet or coil system.

- Phase 1: tethered probe with a magnetic tip; the tether carries power, video and a tool channel.
- Phase 2: untethered capsule, including detecting and freeing a stuck capsule.
- Benchmark: an Endotics-style inchworm probe in the same colon model.
- cfsem: steering field and gradient (force and torque on the tip magnet).
- CortenForge: the full stack: soft colon wall (sim-soft), contact, the tether, controller, meshes.
- Reuses magnet-tracking for localization and wireless-charging to power the capsule.
