# Magnet tracking from sensors (main 1)

Find a permanent magnet's position and orientation from an array of magnetometers.

- cfsem: dipole fields (`cfsem::physics::point_source::dipole`) as the forward model, checked
  against the magnet as current loops and a closed form; coils as the calibration source.
- Solver: the pose from field readings, in plain Rust (CortenForge kept light until its 0.10).
- Hardware (later): cheap 3-axis magnetometer boards; compare the tracked pose with a known path.
  Built simulation-first so the solver and the sensor layout are ready when boards arrive.
- Becomes the localization module for the colonoscopy probe.

Run: `cargo run --release -p magnet-tracking -- <mode>`.

- `check`: the forward models against textbook results (under a second).
- `solve`: the pose solver on a simulated 4 x 4 array: exact recovery, cold starts, noise
  against the Cramer-Rao bound, bias by magnet shape (about a second).

## Checks (`check`, 2026-11-02)

6 of 6 pass. The test magnet: 6 mm across, 6 mm long, Br 1.3 T, tilted.

| Check | Result |
|---|---|
| cfsem's point dipole vs mu0/4pi (3 (m.r) r / r^5 - m / r^3), 200 random points | 1.2e-15 |
| Bulirsch's `cel` vs K(k) and E(k) by the AGM, k 0.1 to 0.9999 | 6.7e-16 |
| Derby and Olbert's closed form on the axis vs the textbook on-axis field, inside and out | 2.3e-10 |
| cfsem current loops vs Derby and Olbert off the axis, 6-60 mm away | 4 loops 2.7e-2, 16 1.8e-3, 64 1.1e-4, 256 7.1e-6 |
| point dipole error, least at L/D 0.87 (theory sqrt(3)/2 = 0.866) | 1.5e-4 at 30 mm |

- Derby and Olbert (Am. J. Phys. 78, 229 (2010)) give a uniformly magnetized cylinder's field
  anywhere through Bulirsch's generalized complete elliptic integral; it is the independent
  reference for the loops (which converge on it at second order, 16x per 4x loops).
- A point dipole is off by 1.8% at 2 diameters from a magnet as long as it is wide (L/D = 1),
  4.3e-3 at 4, 1.1e-3 at 8 and 2.7e-4 at 16 diameters: 4x per doubling, the octupole term. At
  L/D = sqrt(3)/2 that term vanishes; at fixed volume and 30 mm the error falls from 2.8e-3
  (L/D 1.0) to 1.5e-4 (L/D 0.87), against 1.2e-2 at L/D 0.4 and 1.5e-2 at 1.6.
- So a tracking magnet shaped L/D of about 0.87 lets the solver use a point dipole with little
  model error; any other shape needs the loops (or the closed form) near the sensors.

## The pose solver (`solve`, 2026-11-02)

A 4 x 4 array of 3-axis magnetometers at 30 mm pitch (90 x 90 mm, in the plane z = 0); the
magnet 6 mm across and 5.2 mm long (L/D 0.867), Br 1.3 T, moment 0.1521 A m^2; a background like
the Earth's, (20, -5, -40) uT. Nine unknowns from 48 readings: position, axis (two angles),
moment and the background. Levenberg-Marquardt (Marquardt scaling, central-difference Jacobian)
on cfsem's point dipole, in plain Rust (nalgebra). Noise 1 uT per axis is a placeholder until a
real sensor is chosen (stage 3); every error below scales with it.

- No noise, 50 random poses (x and y within 40 mm of centre, 20-80 mm up, any axis), started
  10 mm and about 30 deg off with the moment 10% high and no background: recovered to 3.1e-17 m
  and 0 rad, 11.2 iterations on average.
- Cold start (no prior pose): the best of 162 starts (27 positions x 6 axes) finds all 50 to
  1 um and 1 urad, 7 ms an acquisition. The first try found 34: the other 16 had the axis
  flipped and the moment negative, which gives the same field. The solver now reports the
  positive-moment branch.

With noise, 1000 trials a pose, the spread of the estimates against the Cramer-Rao bound
sigma^2 (J^T J)^-1:

| Pose | x | y | z | axis | bias | z bound, moment known |
|---|---|---|---|---|---|---|
| centre, 30 mm up, axis up | 43.0 / 42.7 um | 42.5 / 42.7 um | 36.0 / 35.7 um | 0.257 / 0.253 deg | 2.0 um | 14.3 um |
| centre, 70 mm up, axis up | 996 / 993 um | 992 / 993 um | 727 / 711 um | 3.30 / 3.30 deg | 24 um | 158 um |
| off the corner, 40 mm up, tilted | 112 / 116 um | 102 / 102 um | 128 / 131 um | 0.396 / 0.401 deg | 3.3 um | 52.5 um |

(spread / bound). The solver reaches the bound everywhere, within the Monte Carlo's own
uncertainty (about 2% for 1000 trials; the largest gap is 3.4%), with a bias far under the
spread.

- Accuracy falls fast with height: from 30 to 70 mm the error grows 23x (about the 3.7th power
  of the distance).
- Knowing the moment helps most in z: the moment and the height trade off (a stronger magnet
  further away reads nearly the same), so fitting the moment costs 2.5x in z at 30 mm and 4.5x
  at 70 mm. A one-time calibration of each magnet's moment is worth it.

A real magnet read by the dipole solver (readings from Derby and Olbert's closed form, no noise):
the position bias, and that bias over the noise's position spread at 1 uT, magnets of the same
volume and Br:

| L/D | 15 mm up | 25 mm | 40 mm | 80 mm |
|---|---|---|---|---|
| 0.5 | 227 um (17x) | 203 um (5.2x) | 155 um (0.96x) | 103 um (0.04x) |
| 0.867 | 5.7 um (0.43x) | 1.6 um (0.04x) | 0.7 um | 0.1 um |
| 1.0 | 75 um (5.6x) | 63 um (1.6x) | 48 um (0.30x) | 32 um (0.01x) |
| 2.0 | 576 um (43x) | 517 um (13x) | 398 um (2.5x) | 265 um (0.11x) |

- With the L/D 0.867 magnet the dipole model's bias stays under half the noise even 15 mm from
  the board; a cube-like magnet (L/D 1) is biased 5.6x the noise there, a long one 43x. Close
  in, any other shape needs the full model (loops or the closed form) in the solver.

## Plan

1. Forward models (done).
2. A simulated sensor array and a pose solver (done): exact recovery, cold starts, the
   Cramer-Rao bound reached, model bias by magnet shape.
3. Accuracy over the workspace for sensor layouts (count, spacing), with a realistic sensor's
   noise and the Earth's field.
4. Tracking a moving magnet along a path (warm starts, speed).
