# Magnet tracking from sensors (main 1)

Find a permanent magnet's position and orientation from an array of magnetometers.

- cfsem: dipole fields (`cfsem::physics::point_source::dipole`) as the forward model; coils as
  the calibration source.
- CortenForge: sensor board and fixture geometry, simulation of the moving magnet, the solver
  (sim-opt) that inverts field readings into a pose.
- Hardware: cheap 3-axis magnetometer boards; compare tracked pose with a known path.
- Becomes the localization module for the colonoscopy probe.
