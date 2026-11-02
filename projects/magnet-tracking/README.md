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
- `layout`: real magnetometers (MMC5983MA, LIS3MDL) with saturation; the weighted solver against
  its bound, 32 layouts over the workspace, and calibration errors (about 3 s).
- `calib`: calibrating a board in simulation: exact recovery, then 40 random boards calibrated
  with jigs of several tolerances, per sensor and by a joint solve, and tracked (about 30 s).

## Checks (`check`, 2026-11-02)

5 of 5 pass. The test magnet: 6 mm across, 6 mm long, Br 1.3 T, tilted.

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

## Sensor layouts (`layout`, 2026-11-02)

Two magnetometers, numbers from their datasheets:

| Sensor | Noise per axis (rms) | Range | Setting |
|---|---|---|---|
| MEMSIC MMC5983MA (Rev. A, 2019) | 0.04 uT (0.4 mG) | +-800 uT (+-8 G) | BW = 00, 50 Hz output |
| ST LIS3MDL (DocID024204 Rev 6) | 0.32 uT x, y; 0.41 uT z (3.2 / 4.1 mG) | +-1200 uT (+-12 G) | ultra-high performance |

- MEMSIC gives 0.4 mG as "total RMS noise"; it is taken per axis here, the pessimistic reading.
  Faster MMC5983MA settings are noisier (0.6 mG and up, to 1000 Hz output). Melexis' MLX90393
  (5-50 mT range) gives its noise only as a plot, so it is left out.
- A sensor reading past its range on any axis is dropped from the fit; each reading is weighted
  by 1 / its noise (the LIS3MDL's z is noisier).
- Check: the weighted solver against the weighted Cramer-Rao bound, 4000 trials, 4 x 4 at
  30 mm. MMC5983MA with the magnet 18 mm up and one sensor saturated (dropped): x 1.07 / 1.07,
  y 1.07 / 1.07, z 0.99 / 0.99 um; LIS3MDL 50 mm up: x 61.0 / 59.9, y 61.7 / 61.4, z 77.7 /
  77.7 um (spread / bound). The first run at 1000 trials read 3-4% over for the LIS3MDL; at 4000
  the largest gap is 1.9%, inside the Monte Carlo's own 1.1%-per-axis uncertainty at under 2 sigma.

The position error bound (moment and background fitted) over a workspace: the magnet anywhere in
x, y within 40 mm of the centre (20 mm steps), axis up, sideways or diagonal; the 95th
percentile at each height:

| MMC5983MA | Board | 20 mm | 40 mm | 60 mm | 80 mm | 100 mm | 120 mm | Most saturated |
|---|---|---|---|---|---|---|---|---|
| 3 x 3 at 30 mm | 60 mm | 28 um | 28 um | 96 um | 286 um | 759 um | 1.8 mm | 2 of 9 |
| 4 x 4 at 30 mm | 90 mm | 4 um | 9 um | 41 um | 136 um | 365 um | 856 um | 2 of 16 |
| 4 x 4 at 40 mm | 120 mm | 7 um | 9 um | 36 um | 110 um | 285 um | 653 um | 2 of 16 |
| 5 x 5 at 30 mm | 120 mm | 3 um | 7 um | 28 um | 89 um | 234 um | 540 um | 2 of 25 |
| 6 x 6 at 30 mm | 150 mm | 2 um | 6 um | 24 um | 71 um | 178 um | 395 um | 2 of 36 |

| LIS3MDL | Board | 20 mm | 40 mm | 60 mm | 80 mm | 100 mm | 120 mm |
|---|---|---|---|---|---|---|---|
| 4 x 4 at 30 mm | 90 mm | 38 um | 77 um | 370 um | 1.2 mm | 3.2 mm | 7.8 mm |
| 4 x 4 at 40 mm | 120 mm | 33 um | 80 um | 321 um | 995 um | 2.6 mm | 6.0 mm |
| 6 x 6 at 30 mm | 150 mm | 23 um | 53 um | 214 um | 643 um | 1.6 mm | 3.6 mm |

(all 32 layouts, 3 x 3 to 6 x 6 at 20-50 mm pitch, in the run's output.)

- The MMC5983MA's 8x lower noise is worth about 8x in error: with it 4 x 4 at 30 mm pitch or
  wider, and every 5 x 5 and 6 x 6, stay under 1 mm (95th percentile) to 120 mm; the LIS3MDL
  reaches 60-80 mm.
- Saturation near the board costs little: at most 5 of 9 sensors (3 x 3 at 20 mm) and usually 2
  drop out; the rest still fix the pose. Small pitches suffer most close in.
- Past 4 x 4 at 30-40 mm the gains are modest (6 x 6 roughly halves the error at 120 mm for 2.25x
  the sensors). Chosen for a first board: 4 x 4 MMC5983MA at 30 mm (90 mm across).

What the noise floor leaves out: calibration. The same board (MMC5983MA, 4 x 4 at 30 mm) with
random gain errors per axis or sensors placed off their nominal positions, no noise, 300 random
boards; position bias, median / 95th percentile (median over the noise's bound):

| Board error | 20 mm up | 60 mm | 120 mm |
|---|---|---|---|
| gain 0.1% per axis | 13 / 27 um (8.6x) | 28 / 50 um (0.9x) | 354 / 596 um (0.5x) |
| gain 1% per axis | 129 / 270 um (86x) | 277 / 502 um (9.0x) | 3.6 / 6.3 mm (4.9x) |
| sensors 0.1 mm off | 158 / 277 um (105x) | 207 / 364 um (6.7x) | 367 / 692 um (0.5x) |
| sensors 0.5 mm off | 791 / 1383 um (526x) | 1.0 / 1.8 mm (34x) | 1.8 / 3.5 mm (2.5x) |

- Calibration, not noise, sets the real accuracy: 1% gain errors or sensors 0.1 mm out of place
  bias the pose 100x the noise floor near the board. Getting to the noise floor needs gains to
  about 0.1% and positions to well under 0.1 mm, so the board must be calibrated against a known
  source (a magnet at known poses, or a coil).
- Far out, gain errors bite through the Earth's field: 0.1% of its 45 uT is a noise-sized false
  signal against the magnet's few uT, so the bias grows again at 120 mm.

## Calibration (`calib`, 2026-11-02)

Each sensor reads g (B + E) + o: a 3 x 3 gain matrix g (gain and small axis misalignment), an
offset o, the background E, at a position d off its nominal one; 15 unknowns a sensor (g, c =
g E + o, d), 240 for the board. A jig sets the calibration magnet (the same 6 x 5.2 mm one) at
196 poses: a 7 x 7 grid at 15 mm over the board, 25 and 45 mm up, axis up or along x. Readings
come from the exact cylinder (Derby and Olbert), so the magnet's shape error is in them; the
fit models it as a dipole. Readings are averaged 16x (noise / 4); saturated ones are skipped.

- Two degeneracies shape the fit. A gain common to every sensor looks exactly like a different
  magnet moment, which the tracker fits anyway, so the calibration magnet's moment is held at its
  nominal value. The background at calibration folds into c; the tracker's fitted background
  takes any change after.
- Check: no noise, a perfect jig, dipole readings, a board with 1% gains, 1 deg misalignment,
  0.5 uT offsets and 0.5 mm position errors: every sensor's unknowns recovered to 2.2e-16 (gain),
  1.4e-20 T (offset), 4.9e-18 m (position).

40 random boards with those errors (gains 1%, misalignment 1 deg, offsets 0.5 uT, positions
0.5 mm: assumed for a hand-built board), calibrated, then tracking 9 test poses without noise;
position bias, median / 95th percentile, and the noise floor's bound:

| Calibration | 20 mm up | 60 mm | 120 mm |
|---|---|---|---|
| none | 825 um / 1.56 mm | 1.47 / 2.51 mm | 13.0 / 24.3 mm |
| perfect jig | 1.4 / 1.8 um | 1.2 / 3.1 um | 27.7 / 51.7 um |
| jig off 0.05 mm (position only) | 33.8 / 66.2 um | 99.1 / 168 um | 0.92 / 1.79 mm |
| jig off 0.2 deg (angle only) | 23.0 / 44.0 um | 65.7 / 122 um | 1.01 / 1.91 mm |
| jig off 0.05 mm and 0.2 deg | 43.5 / 77.1 um | 127 / 201 um | 1.29 / 2.89 mm |
| jig off 0.2 mm and 0.5 deg | 161 / 267 um | 461 / 762 um | 4.53 / 8.67 mm |
| same jig, joint solve | 31.7 / 63.6 um | 32.3 / 61.8 um | 69.6 / 108 um |
| jig off 0.5 mm and 1 deg, joint solve | 73.9 / 144 um | 83.8 / 153 um | 157 / 229 um |
| (noise floor) | 2.1 um | 34.2 um | 790 um |

- Calibration against a perfect jig takes an uncalibrated board's 0.8-13 mm down to under the
  noise floor everywhere: the 240 unknowns are all recoverable from 196 poses.
- Taking the jig's poses as exact, its errors go straight into the board: 0.05 mm of position or
  0.2 deg of angle error each cost tens of um near the board and about 1 mm at 120 mm (they add
  in quadrature). A printed jig (0.2 mm, 0.5 deg) would leave 0.16-4.5 mm.
- The fix: solve every pose and every sensor at once (bundle adjustment: 196 x 5 + 16 x 15 =
  1,220 unknowns), with the jig's nominal poses as measurements of known tolerance rather than
  exact, and the sensors' nominal positions with 1 mm. The readings alone leave the whole
  arrangement free to move, turn and scale (a scaling's c^-3 in the fields is a common gain); the
  jig's priors fix that with their errors averaged over 196 poses. With the 0.2 mm, 0.5 deg jig
  the bias falls 5-65x, to 32 um at 20 and 60 mm and 70 um at 120 mm, under the noise floor from
  60 mm out; even a 0.5 mm, 1 deg jig gives 74-157 um.
- What is left is mostly where the whole board sits in the jig's frame: its mean error over 196
  poses, about 0.2 mm / sqrt(196) = 14 um an axis (25 um in 3D) for the 0.2 mm jig, 36 um an axis
  for the 0.5 mm one. Positions relative to the board, and any motion, don't see it.
- Getting there took three tries. Alternating (sensors given poses, then poses given sensors,
  re-anchored to the jig each round) with a rigid anchor gave 58 / 173 / 1400 um; adding scale to
  the anchor, 48 / 90 / 665 um at 20 rounds and 38 / 49 / 246 um at 60, still falling: slow
  convergence along the coupled directions. The joint solve converges in a few dozen steps
  (the whole run about 30 s).

## Plan

1. Forward models (done).
2. A simulated sensor array and a pose solver (done): exact recovery, cold starts, the
   Cramer-Rao bound reached, model bias by magnet shape.
3. Sensor layouts (done): two real magnetometers, saturation, 32 layouts; a 4 x 4 MMC5983MA
   board at 30 mm stays under 1 mm to 120 mm; calibration, not noise, limits a real board.
4. Calibration (done): 15 unknowns a sensor from a magnet in a jig; with a perfect jig, under the
   noise floor; a joint solve (bundle adjustment) makes a printed jig good to 32-70 um.
5. Tracking a moving magnet along a path (warm starts, speed).
