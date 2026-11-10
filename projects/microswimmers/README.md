# Helical microswimmers (extra)

Tiny corkscrews that swim when a magnetic field turns them. The long-term hope behind robots
like these is chemotherapy that goes only to the tumour instead of through the whole body. This
project shows the physics on a desk, in the open, with every number from a real run, and says
plainly what is still unsolved.

## Why this matters

Chemotherapy travels in the blood to every organ, and that is where most of its side effects
come from. Researchers are building swimmers too small to see without a microscope that carry a
drug, are steered from outside the body by magnetic fields, and release it at the tumour.

Where it stands, honestly:

- Shown in animals, not approved in people. In a well-known study, magnetic bacteria carried
  drug-filled capsules into the oxygen-starved cores of tumours in mice, where drugs struggle to
  reach (Felfoul et al., Nature Nanotechnology 11, 941 (2016)).
- Nearby magnetic medicine is already in clinics: magnetically steered heart catheters
  (Stereotaxis, since the early 2000s), a magnetically steered capsule camera for the stomach
  (NaviCam, FDA-cleared in 2020) and magnetic nanoparticles that heat brain tumours (NanoTherm,
  approved in the EU in 2010).
- Lab swimmers like ours date from 2009: artificial bacterial flagella (Zhang et al., Appl.
  Phys. Lett. 94, 064107 (2009)) and nanostructured propellers (Ghosh and Fischer, Nano Lett.
  9, 2243 (2009)).

What this project is not: a medical device. No biology, no claims about treating anyone. It is
the physics all of the above depends on, made visible.

## Why corkscrews, in plain words

At the size of a cell, water acts like thick honey. Inertia is gone: the moment a swimmer stops
pushing, it stops. A back-and-forth stroke, like a paddle, undoes itself on the way back
(Purcell's scallop theorem, "Life at low Reynolds number", Am. J. Phys. 45, 3 (1977)). Bacteria
get round it by spinning a corkscrew tail. Spinning only ever goes one way, and the helix turns
it into forward motion, like a screw into wood.

A magnetic swimmer copies that. A small magnet in the swimmer, magnetized across its axis,
follows a field that rotates. The swimmer turns and the helix screws it forward. Tilt the plane
the field rotates in and the swimmer turns with it: that is steering.

Why twist instead of pull: pulling a magnet needs a field that changes over distance (a
gradient), and gradients from coils outside the body fade fast as the workspace grows. A
rotating, uniform field only has to twist the swimmer, which scales much better to a
body-sized workspace (Abbott et al., "How should microrobots swim?", Int. J. Robotics Research
28, 1434 (2009)).

## The desk demo

Cell-sized swimmers can't be made at home, and they don't need to be. What decides the physics
is the ratio of inertia to stickiness, the Reynolds number. A printed corkscrew a couple of
centimetres long in glycerol is in nearly the same regime as a microrobot in water: both are
under 1, the desk demo at 0.05-0.2 and a cell-sized swimmer at about 0.001 (`swim`, below).
So:

- a printed helix a few mm long with a tiny magnet magnetized across its axis,
- a jar of glycerol inside three pairs of coils (x, y, z),
- two pairs driven 90 degrees apart make the rotating field; the third tilts it to steer,
- a camera above, to measure speed against the model and to steer through a printed maze.

Something people can watch and understand is easier to trust than a promise.

## Plan

`cargo run --release -p microswimmers -- <mode>`, modes added as they are built (so far:
`check`, `swim`, `coils`, `shape`, `steer`).

1. `check` (done, results below): the physics before anything is built on it.
   - The reference model: regularized Stokeslets (Cortez, Fauci and Medovikov, Phys. Fluids 17,
     031504 (2005)), points on the body's surface, in plain Rust. Checked against Stokes' sphere
     and the exact prolate spheroid (Chwang and Wu, J. Fluid Mech. 67, 787 (1975)), whose
     formula is itself checked against the sphere and slender-body theory.
   - The quick model: resistive force theory (Gray and Hancock, J. Exp. Biol. 32, 802 (1955);
     Lighthill, SIAM Rev. 18, 161 (1976)), mapped against the reference over pitch angle and
     wire thickness.
   - Chirality: a straight wire spun about its axis doesn't swim; a mirrored helix swims
     backwards.
   - CortenForge's own fluid model (`cortenforge-sim-core` 0.9, `src/forward/passive.rs`): its
     viscous drag is the same in every direction, so a helix built from its geoms can't swim in
     syrup. Confirmed by a run and logged in docs/CORTENFORGE.md with a Suggestion. The
     swimmer's fluid forces come from our model, applied through `xfrc_applied`.
   - cfsem: a Helmholtz pair's centre field, and the field of two pairs in quadrature.
2. `swim` (done, results below): speed against rotation rate for a printable swimmer in
   glycerol, up to the step-out frequency, where the magnet can't twist it fast enough and it
   slips and slows. Whether it sinks faster than it swims. The Reynolds numbers of the desk demo
   and a cell-sized swimmer in water side by side.
3. `coils` (done, results below): three Helmholtz pairs round the jar (coil-designer's coil
   model): whether they fit, how even the rotating field is over the jar, inductance, the drive
   voltage at the rotation rates `swim` needs, power and heating, and the leftover gradient's
   pull against the twist.
4. `shape` (done, results below): the helix (radius, pitch angle, wire, turns) chosen for the
   fastest climb at 2 mT with sim-opt (RicherSA, a grid for reference), printable in resin; then
   the part as a cf-design solid, meshed, checked with mesh-printability (SLA) and written as STL.
5. `steer` (done, results below): sim-core with gravity, buoyancy and contact with the jar and a
   printed maze; a controller tilts the rotation plane to follow a path (sim-rl if a hand-tuned
   one isn't enough; it was enough).
6. `scale`: the honest numbers. Coils sized for a human torso, their power and stored energy,
   and the swimmer's speed against blood flow. This fills the section below with real figures.
7. Hardware: print the swimmer and the coil formers, a three-channel driver, the jar and a
   camera; measure speed against rotation rate and compare with `swim`.

Connects to coil-designer (coils and formers), magnet-tracking (open question: can
magnetometers find the swimmer's magnet if the drive field pauses while they read?) and
colon-probe (the same steer-from-outside idea at a larger size).

## Checks (`check`, 2026-11-09)

24 of 24 pass in 1-2 minutes. In plain words: our model of a swimmer in syrup reproduces exact
textbook results to about 1%. The quick textbook shortcut gets either the speed or the twisting
effort right, never both. CortenForge's built-in fluid model can't make a corkscrew swim at all,
so the swimmer's forces come from our model.

**The regularized Stokeslet** matches the plain Stokeslet far away (5.0e-5 at 100 eps, 5.0e-7
at 1000, falling as (eps/r)^2). The flow is divergence-free to 2.2e-9 (central differences),
and the velocity at the blob's centre is exactly f / (4 pi mu eps).

**A sphere**, eps = c h, h the spacing between points. Errors of drag / torque, % (three of the
five spacings run):

| h / a | points | c = 0.1 | c = 0.2 | c = 0.3 | c = 0.4 | c = 0.5 | c = 0.75 |
|---|---|---|---|---|---|---|---|
| 0.785 | 22 | -38.8 / -54.9 | -16.4 / -25.6 | -4.5 / -4.2 | +3.3 / +13.1 | +9.2 / +28.1 | +21.2 / +63.5 |
| 0.393 | 82 | -25.5 / -39.7 | -9.8 / -15.9 | -2.9 / -2.8 | +1.2 / +6.0 | +4.1 / +12.7 | +9.4 / +26.2 |
| 0.100 | 1242 | -8.2 / -14.6 | -2.8 / -4.7 | -0.84 / -0.77 | +0.20 / +1.46 | +0.87 / +3.0 | +2.0 / +5.7 |

c = 0.3 has the least worst error at every spacing and is used from here on (`EPS_C`, which
`check` confirms). (h = 0.785 a is a
tube with 8 points round the wire.)

**Prolate spheroids.** The exact formulas agree with the sphere to 1.2e-6 (eccentricity 1e-3)
and with slender-body theory's leading terms to 6.6e-6 (aspect 1000). The model's errors, %
(force along / force across / torque about the axis):

| aspect | 8 round | 12 round | 16 round | extrapolated |
|---|---|---|---|---|
| 2 | -3.2 / -5.4 / -7.2 | -1.9 / -3.6 / -4.2 | -1.5 / -2.8 / -3.0 | -0.33 / -0.62 / +0.45 |
| 5 | -1.2 / -4.0 / -5.6 | -1.1 / -3.0 / -3.5 | -0.73 / -2.3 / -2.3 | +0.20 / +0.08 / +1.25 |
| 10 | -0.65 / -3.3 / -4.7 | -0.66 / -2.5 / -3.2 | -0.49 / -1.9 / -2.3 | +0.03 / -0.15 / +0.27 |

- The error is always low and falls about as the spacing. 16 points round alone is up to 3.02%
  low: the first run's check (3%) failed on it.
- Extrapolating from 12 and 16 round (R = 4 R16 - 3 R12, for an error proportional to h) brings
  the worst to 1.25%. The extrapolated model is the reference from here on.

**Chirality.** A mirrored helix advances backwards (sum 5.5e-12 of the advance), a straight wire
spun about its axis makes no thrust (3.7e-14 of the helix's), and the resistance matrix is
symmetric to 2.2e-12. A right-handed helix turned counter-clockwise (seen from +z) moves to +z,
like a screw.

**A helix converging** (45 deg, wire b = 0.2 R, 2 turns): advance per turn as a fraction of the
pitch 0.1916 / 0.1968 / 0.1993 at 8 / 12 / 16 points round, extrapolated 0.2067. The torque,
T / (mu omega R^2 L), goes 4.197 / 4.236 / 4.256, extrapolated 4.315. Two turns wobble: the axis
tilts at 5.3% of the spin rate.

**Resistive force theory against the model** (extrapolated), 2 turns, ratios:

| pitch angle | b / R | points (16 round) | advance / pitch | Gray-Hancock | Lighthill | torque | Gray-Hancock | Lighthill |
|---|---|---|---|---|---|---|---|---|
| 30 | 0.1 | 10322 | 0.1364 | 0.968 | 1.048 | 3.023 | 0.597 | 0.939 |
| 30 | 0.2 | 5202 | 0.1244 | 0.997 | 1.059 | 4.189 | 0.557 | 0.911 |
| 45 | 0.1 | 7330 | 0.2317 | 0.965 | 1.045 | 3.007 | 0.581 | 0.912 |
| 45 | 0.2 | 3698 | 0.2067 | 1.005 | 1.064 | 4.315 | 0.516 | 0.869 |
| 60 | 0.1 | 6002 | 0.2895 | 0.993 | 1.100 | 2.894 | 0.609 | 0.901 |
| 60 | 0.2 | 3042 | 0.2468 | 1.061 | 1.161 | 4.259 | 0.522 | 0.849 |

- Gray and Hancock's coefficients get the speed within 6.1% but the torque 39-48% low.
  Lighthill's get the torque 6-15% low and the speed 4.5-16% high.
- Both drift most for thick wire at steep angles, where printed swimmers live. So `swim` uses
  the model, and resistive force theory is a quick estimate only.

**CortenForge's fluid model** (sim-core 0.9; glycerol, 1.41 Pa s, fluid density 0):

- A free 1 mm sphere with `fluidshape="ellipsoid"` settles to Stokes' drag and torque exactly
  (1.6e-14, 4.7e-15). With the default body model (an inertia box) it is 22.5% low on drag and
  53.5% low on torque, as that model's formula gives for a sphere (its box is 1.549 a a side).
- A sphere pushed through the last three entries of `xfrc_applied` moves without turning, so the
  layout is [torque, force], confirmed by a run.
- The desk helix (R 1.5 mm, pitch 6 mm, 57.5 deg, wire 0.3 mm, 2 turns) in the model advances
  0.2455 of its pitch a turn (1.47 mm). Built in sim-core from 32 capsules and twisted the same
  way, it doesn't advance at all: 5e-18 of the pitch a turn with the ellipsoid model, 2e-21 with
  the inertia box. Its turning drag is 2.61x the model's (ellipsoid) and 3.81x (inertia box):
  each piece is dragged as if alone, with no shielding by its neighbours. Glycerol's density
  (1261 kg/m^3) adds an advance under 2e-7, of either sign: noise from the fast-flow terms.

**cfsem**: the Helmholtz centre field is within 3.3e-12 of (4/5)^1.5 mu0 NI / R. Two pairs in
quadrature (one nested at 1.2 R, its current raised to match) give a centre field that turns at
omega t with constant size (3.3e-16).

## Swimming (`swim`, 2026-11-09)

In plain words: the desk swimmer is about 2.3 cm long and moves a fifth of its pitch each turn,
1.2 mm/s for every turn a second. A 5 mT rotating field can turn it up to 9 times a second, for
10.6 mm/s, nearly half a body length a second. Turn the field faster and the magnet slips and the
swimmer slows down. It's heavier than glycerol and sinks at 3-4 mm/s, so it needs a field of at
least about 2 mT to swim upwards. About 17 s.

**The swimmer** (`Swimmer::desk`):

- A resin helix (1150 kg/m^3): R 1.5 mm, pitch 6 mm (57.5 deg), wire 0.6 mm across, 3 turns.
- A capsule head 2.5 mm across and 4.5 mm long, 0.5 mm ahead of the helix, holding a 2 x 2 mm
  magnet magnetized across its axis. The magnet is N42-class (Br 1.3 T, as in magnet-tracking),
  7500 kg/m^3, m = 6.5e-3 A m^2.
- A thin strut joins head and helix; its drag is ignored.
- In all: 23.5 mm long, 27.6 mm^3, 71.6 mg.
- The model: 7183 points at 16 round, extrapolated from 12 and 16.

**Per turn:**

| | advance / pitch | mm/s per Hz | torque, N m per Hz | wobble |
|---|---|---|---|---|
| axis held by the field | 0.1973 | 1.184 | 3.624e-6 | 0 |
| axis free to wobble | 0.1952 | 1.171 | 3.606e-6 | 0.018 |

- At 12 and 16 round alone the advance is 0.1869 and 0.1895. The extrapolation adds 4%, and the
  extrapolated model was within 1.25% on the spheroids (`check`).
- The real swimmer lies between held and free, because the field also resists tilting. The two
  differ by 1%; held is used below. (`steer` later found the free axis is the right one: sim-core
  matches it to 0.1%.)
- Below step-out the advance per turn depends only on the shape, not on the fluid (Stokes flow is
  linear). Viscosity sets only the step-out rate and the sinking speed. Glycerol's viscosity
  changes a lot with temperature and the water it soaks up from the air, but the measured
  advance per turn shouldn't.

**Step-out**, f = m B / (2 pi D), D the turning drag:

| B, mT | step-out, Hz | top speed, mm/s | body lengths / s |
|---|---|---|---|
| 1 | 1.79 | 2.12 | 0.091 |
| 2 | 3.59 | 4.25 | 0.181 |
| 5 | 8.97 | 10.62 | 0.453 |
| 10 | 17.94 | 21.23 | 0.905 |

**Past step-out** (5 mT, Adler's equation): the magnet slips a turn now and then, and the speed
falls from 10.62 mm/s at step-out to 6.81 at 1.1 times the rate, 4.06 at 1.5x, 2.84 at 2x and
1.82 at 3x. Speed against field rate rises in a straight line, peaks at step-out and drops. The
peak is a clean thing to measure with a camera.

**Gravity:**

- In glycerol it weighs 36.8 mg. It sinks at 3.04 mm/s with its axis held level and 3.82 mm/s
  held upright.
- To outswim sinking, the field must beat 1.43 mT when level. To climb straight up, it must beat
  1.80 mT and the swimmer must turn faster than 3.22 Hz.
- The heavy head puts the centre of mass 2.33 mm ahead of the centre of buoyancy. Gravity's
  tilting torque is 5.0% of what the magnet can hold at 5 mT.
- A lighter head (hollow, or foam round the magnet) would sink slower: a lever for `shape`.

**Reynolds numbers**, by the whole length:

- At step-out, the desk demo is at 0.045 (1 mT), 0.089 (2 mT), 0.22 (5 mT) and 0.45 (10 mT);
  turning, 0.11 at 5 mT.
- The same shape 50 um long in water at 10 Hz (25 um/s) is at 1.3e-3.
- So the desk demo is in the same sticky regime at modest fields. At the 5 mT top speed, a rough
  size of inertia's share is 4.2%: Oseen's first correction, (3/16) Re, for a sphere as wide as
  the swimmer is long, which overestimates it. At 10 mT it is twice that, so keep the demo slow
  or the fluid thick.
- The swimmer's own inertia settles in 0.76 ms, far quicker than a turn.

Not yet in the model: the jar's walls (drag rises near a wall), the strut, and how strongly the
field holds the axis.

## Coils (`coils`, 2026-11-09)

In plain words: three pairs of coils, 15, 18 and 22 cm across, nest round a 56 mm jar. 2 mT is
enough for the swimmer to climb. It takes 7-13 W and up to 11.7 V, so a 12 V supply is marginal
and 24 V comfortable. 5 mT takes 46-84 W and up to 30 V, and would warm the copper by up to
7.6 C a minute with no cooling at all: short bursts only. The field is even: across all of the
glycerol its average strength stays within 1.4% of the centre's (2.5% below at its weakest moment
in a turn), and its axis tilts under 1 degree.
The magnet's pull from the field's unevenness is at most 6% of the swimming speed, and 0.2% in
the middle. About 30 s.

**The coil set.** Each coil is AWG 20 (0.812 mm bare), 13 x 13 = 169 turns in an 11.3 mm square
pack, two coils a pair at Helmholtz spacing (the mean radius). It uses coil-designer's model:
one cfsem filament per turn, inductance and resistance at 20 C, checked there (58 of 58).

| pair | mean radius, mm | mT/A | ohm | mH | L/R, ms | wire, m | copper, kg |
|---|---|---|---|---|---|---|---|
| x | 75.0 | 2.0254 | 5.305 | 17.089 | 3.221 | 159.3 | 0.739 |
| y | 92.0 | 1.6513 | 6.508 | 22.300 | 3.427 | 195.4 | 0.906 |
| z (upright) | 112.0 | 1.3566 | 7.923 | 28.720 | 3.625 | 237.9 | 1.103 |

**Fit.**

- The windings clear each other by 3.9 mm (x to y), 7.2 mm (y to z) and 26.2 mm (x to z).
- The jar (56 mm across, 60 mm tall) stands 35-81 mm from every winding. Lifting it straight out
  passes 3.8 mm from the x windings and 12.3 mm from the y windings, so the x pair can't shrink
  much if the jar is to come out of the top.
- Smaller coils would cost less power (at a fixed pack, power grows as the cube of the size),
  but they would have to be built round the jar.
- Orthogonal pairs on one centre have no mutual inductance (each pair's field through the
  others' loops cancels by symmetry), so the three drive independently.

**Driving it.** The rate is the swimmer's step-out (3.59 Hz at 2 mT, 8.97 Hz at 5 mT). Peak
current and voltage are per pair (x, y, z), with the mean power of all three. Heating is for the
worst pair with no cooling at all (an upper bound).

| B, mT | swimming | peak A | peak V | mean W | worst C/min | peak A/mm^2 |
|---|---|---|---|---|---|---|
| 2 | climbing (field turns in x-y) | 0.99, 1.21, 0 | 5.3, 7.9, 0 | 7.4 | 0.82 | 2.34 |
| 2 | level along x (y-z) | 0, 1.21, 1.47 | 0, 7.9, 11.7 | 13.4 | 1.22 | 2.85 |
| 2 | level along y (x-z) | 0.99, 0, 1.47 | 5.3, 0, 11.7 | 11.2 | 1.22 | 2.85 |
| 5 | climbing | 2.47, 3.03, 0 | 13.3, 20.1, 0 | 46.0 | 5.13 | 5.85 |
| 5 | level along x | 0, 3.03, 3.69 | 0, 20.1, 29.8 | 83.6 | 7.60 | 7.12 |
| 5 | level along y | 2.47, 0, 3.69 | 13.3, 0, 29.8 | 70.0 | 7.60 | 7.12 |

- L/R is 3.2-3.6 ms against a turn of 110 ms or more, so at these rates the coils are nearly pure
  resistors: inductance adds under 2% to the voltage.
- Copper's resistance rises 0.39% per degree, so a warm coil needs more voltage for the same
  field.
- Heat that reaches the jar thins the glycerol. That raises step-out and sinking, but not the
  advance per turn (`swim`).

**The field over the glycerol**, per unit B0, on 1469 points filling the 50 x 50 mm of glycerol
(the core is the middle 10 mm):

| swimming | where | mean size, min-max | weakest | ellipticity | axis tilt, deg |
|---|---|---|---|---|---|
| climbing | all | 0.992-1.013 | 0.975 | 0.036 | 0.94 |
| climbing | core | 1.000-1.000 | 1.000 | 0.000 | 0.01 |
| level along x | all | 0.997-1.007 | 0.989 | 0.015 | 0.35 |
| level along x | core | 1.000-1.000 | 1.000 | 0.000 | 0.00 |
| level along y | all | 0.990-1.014 | 0.975 | 0.031 | 0.79 |
| level along y | core | 1.000-1.000 | 1.000 | 0.000 | 0.01 |

- The weakest point in a turn sets step-out there: it is 2.5% lower at worst.
- Ellipticity is how far the field's path is from a circle (largest over smallest, less 1).
  Tilt is the angle between the field's turning axis and the way the swimmer should go: the
  steering error, under 1 degree anywhere.

**The pull.** The magnet is drawn towards a stronger field. The force is (m . grad) B with the
magnet trailing the field by its lag (sin phi = (f / f_so) B0 / |B|), averaged over a turn. The
drift it causes is that force through the swimmer's drag, the swimmer pointing the way it swims.
At 0.75 of step-out the worst drift is 5.6% / 3.0% / 5.9% of the swimming speed over all the
glycerol (climbing / level along x / level along y), and 0.2% / 0.1% / 0.2% in the core. Force
and swimming speed both scale with B0, so the ratio holds at any field.

Not yet: the coil formers (coil-designer's tooling), the driver electronics, and the field
measured with a Hall probe.

## Shape (`shape`, 2026-11-09)

In plain words: a printable head for the magnet has to be thicker than `swim`'s, and on that
head the first design barely climbs at 2 mT (0.15 mm/s). The search found a longer, more open
corkscrew of thinner wire: two turns, 3.6 mm across, a 12.2 mm pitch, 0.51 mm wire and 28.8 mm
long. It climbs at 2.12 mm/s at 2 mT, needs only 1.27 mT to climb at all, and tops out 50%
faster (5.84 mm/s at 2 mT, 14.6 at 5 mT). The part meshes watertight, and the resin check finds
nothing it shouldn't. The search takes about 3 min; the STL, about 15.

Modes: `shape` (the search, writes `out/shape/evals.csv` and `best.txt`); `shape part [cell mm]
[octree [snap steps]]` (the part for `best.txt`, written to `out/shape/swimmer.stl`); `shape
probe` (where the octree mesher's defects come from, below). `out/` is git-ignored.

**The head.** CortenForge's SLA preset (`PrinterConfig::sla_default`) wants walls of at least
0.4 mm, so the 2 x 2 mm magnet (0.05 mm clearance) sits in a head 2.9 mm across. `swim`'s head
was 2.5 mm, with 0.25 mm walls. The magnet goes in through the nose, so the nose is flat and
open (sealed after); the model now has that flat nose too.

**The search.** The score is the climb at 2 mT: top speed (at step-out) less sinking, axis
upright.

- Limits: neighbouring turns at least 0.5 mm apart, at most 30 mm long.
- Box: helix radius 1-3 mm, pitch angle 35-70 deg, wire 0.5-1.0 mm across; 2, 3 or 4 turns.
- Each evaluation uses the model at 6 and 8 points round the wire, extrapolated (0.7 s). On the
  baseline it reads 1.8% low on top speed and equal on sinking against the full model (12 and
  16 round).
- CortenForge's RicherSA (sim-opt, 60 evaluations, one seed) runs for each number of turns,
  with a 4 x 4 x 3 grid as a reference.

| turns | RicherSA best, climb mm/s | grid best, climb mm/s | the better, full model |
|---|---|---|---|
| 2 | +2.10 (R 1.56 mm, 38.9 deg, wire 0.51 mm) | +1.86 (R 1.67 mm, 46.7 deg, 0.50 mm) | +2.12 |
| 3 | +0.76 (R 1.42 mm, 50.7 deg, 0.64 mm) | +0.66 (R 1.00 mm, 46.7 deg, 0.50 mm) | +0.84 |
| 4 | -0.74 (R 2.10 mm, 66.8 deg, 0.51 mm) | +0.96 (R 1.00 mm, 46.7 deg, 0.50 mm) | +1.02 |

- RicherSA beat the grid at 2 and 3 turns but missed badly at 4, stuck at a wide, steep helix
  that steps out at only 0.84 Hz per mT. Keeping the better of the two guards against that.
  With one seed and 60 evaluations, the winner is the best found, not a proven optimum.
- The winner sits near the limits: the thinnest wire allowed, nearly the longest swimmer, and a
  low pitch angle. Printing thinner wire or allowing a longer swimmer would buy more.

**The chosen swimmer** (full model) against the baseline (`swim`'s helix on the printable head):

| | baseline | chosen |
|---|---|---|
| helix | R 1.5 mm, pitch 6 mm, 0.6 mm wire, 3 turns | R 1.56 mm, pitch 12.19 mm, 0.51 mm wire, 2 turns |
| length, mass | 22.5 mm, 73.8 mg | 28.8 mm, 70.2 mg |
| advance a turn / pitch | 0.1909 | 0.1358 |
| step-out | 1.696 Hz per mT | 1.764 Hz per mT |
| top speed at 2 mT | 3.89 mm/s | 5.84 mm/s |
| sinking upright / level | 3.74 / 3.03 mm/s | 3.71 / 2.85 mm/s |
| climb at 2 mT | +0.15 mm/s | +2.12 mm/s |
| least field to climb | 1.92 mT | 1.27 mT |

It advances less of its pitch a turn, but its pitch is twice as long, so it covers more ground a
turn for about the same twisting effort.

**The part** (cf-design, mm): the helix swept with `pipe_spline`, a strut of the same wire from
its end to the head, the head as a capsule cut flat, and the magnet pocket open at the nose.

- Meshed with the uniform marching-cubes mesher at 0.051 mm cells, welded at 10 nm: 114,156
  faces, watertight and manifold, in 915 s.
- Volume 19.39 mm^3, against the model's 19.23 mm^3 less the pocket; the strut adds the rest.
- SLA check, printed upright: every thin-wall flag is on the head's wall, at 0.397-0.398 mm (the
  0.4 mm wall as designed, less faceting), none critical.
- Two overhang regions are critical, the largest 10.3 mm^2 at up to 83 deg. Resin printing
  handles overhangs with supports, which this check doesn't model.
- The check also flags 92 overlapping triangle pairs. On a design differing only in the fifth
  digit (`shape probe`) it flagged 31, so these look like the check's known touching-sliver
  false positives (docs/CORTENFORGE.md); not inspected.

**Why the uniform mesher** (`shape probe`, 0.051 mm cells). CortenForge's fast octree mesher
folds triangles on the thin curved sweep:

| piece | octree: overlap pairs, thin walls under 0.2 mm | marching cubes: same | time, octree / marching cubes |
|---|---|---|---|
| capsule, 10 mm, wire radius | 0, 8 (0 snapped) | 0, 0 | 0.0 / 0.0 s |
| straight `pipe_spline`, 10 mm | 0, 8 (0 snapped) | 0, 0 | 1.5 / 9.3 s |
| the helix wire | 101, 4681 (0.000 mm) | 26, 0 | 6.0 / 695 s |
| the whole part | 102, 4815 (0.000 mm) | 31, 0 | 155 / 918 s |

Snapping vertices onto the surface doesn't fix the helix, so the folds come from the mesher.
Logged in docs/CORTENFORGE.md (`cf-design`).

Not yet: a print. The strut is in the part but not in the model, and the glue that seals the
nose is in neither.

## Steer (`steer`, 2026-11-09)

Why this step matters: swimming isn't enough. The point of a drug-carrying swimmer is to take the
medicine to one place. Chemotherapy is hard on people because the drug reaches fast-dividing
cells all over the body, not just the tumour's: the lining of the gut, hair, the bone marrow that
makes blood. A swimmer that can be steered to a tumour could release the drug there instead.
The printed staircase below stands in for that trip: a winding route to a target, with walls to
miss.

What it leaves out: here a camera sees the swimmer through clear glycerol, which a body doesn't
allow. Finding the swimmer inside a person stays on the unsolved list (below).

In plain words: in CortenForge's physics engine, the swimmer `shape` chose stands itself up off
the jar's floor, laps a printed post and threads two windows 20 mm apart in height, in 21 s at
5 mT. Steering takes more field than climbing. At 2 mT it climbs straight up but can't pitch up
from level; the course needs 4 mT (38 W) or, better, 5 mT (62 W, the worst coil warming 1.5 C
over a run). On the way: the heavy magnet in the nose makes it sink nose first like a dart, it
flies tilted against its own sinking like a helicopter, a simple controller beat a cleverer one,
and two contact bugs turned up in sim-core. About 1 min.

Modes: `steer`, `steer course` (the course alone), and `steer wall` and `steer contact`
(sim-core's contacts, below).

**The swimmer in sim-core.**

- A free body with the model's mass (70.2 mg) and centre of mass. Its inertia comes from samples
  (their mass within 0.1% of the model's), and gravity is sim-core's own.
- Buoyancy acts at the centre of buoyancy. The magnet adds its twist, m x B, and its pull,
  (m . grad) B. The fluid's force comes from our 6 x 6 resistance matrix. All three go in
  through `xfrc_applied`.
- The swimmer settles to a new speed in 0.55 ms and a new spin in 0.10 ms, so a plain drag force
  would need steps of microseconds. Instead each step applies the drag of the velocity the step
  will end with (linearly implicit). Sim-core lands on that velocity to 2e-4 at 1 ms steps, and
  its mass matrix matches ours to 1e-20.
- For contact the swimmer is 143 spheres along the wire, the strut and the head (why spheres:
  below).

**Checked against the model** (uniform 2 mT field, 0.75 of step-out, 2.68 Hz):

| | sim-core | model, axis free | model, axis held |
|---|---|---|---|
| swimming alone (weightless), mm/s | 4.336 | 4.331 | 4.438 |
| sinking alone, no field, mm/s | 3.760 | 3.759 (free to turn) | 3.715 |
| climbing (both), mm/s | 0.619 | 0.617 | 0.724 |
| spin past step-out at 1.5x (Adler), Hz | 1.364 | 1.366 | |

- The swimming speed changes by 0.06% from 0.02 ms steps to 4 ms steps.
- The axis is free to wobble. The field resists slow tilting but not the wobble at each turn.
  `swim` and `shape` used the held axis, which is 2.5% optimistic on speed for this swimmer (1%
  for `swim`'s). `shape`'s choice stands.
- Released level with no field, it turns nose down in 11.3 s and then sinks spinning slowly
  (0.026 Hz): a corkscrew turns as it falls.

**Why steering takes more field than climbing.**

- The magnet puts the centre of mass 2.8 mm from the nose and 26 mm from the tail, and drag on
  sideways motion acts 10 mm behind it. Sinking pulls the nose down, like a dart.
- Swimming level, the nose sits below the field's axis and swings to the side. The turning field
  turns part of a tilting torque into a sideways one. Settled for 30 s:

| B, mT | rate, Hz | nose down, deg | swung sideways, deg | in all, deg | forward, mm/s | sinking, mm/s |
|---|---|---|---|---|---|---|
| 2 | 2.68 | 21.2 | 31.6 | 37.4 | 3.70 | 4.53 |
| 3 | 4.02 | 15.5 | 20.0 | 25.1 | 6.11 | 4.64 |
| 4 | 5.36 | 11.9 | 14.8 | 18.9 | 8.38 | 4.67 |
| 5 | 6.70 | 9.7 | 11.8 | 15.2 | 10.62 | 4.69 |

- It sinks at 4.5-4.7 mm/s rather than the 2.85 of a swimmer held level, because the drooping
  nose also swims it downwards.
- Pitching up from level, with the field's axis swung up at a steady rate:
  - at 2 mT it loses hold and tumbles even at 0.05 rad/s;
  - 3 mT keeps up to 0.1 rad/s, 4 mT to 0.3, and 5 mT to 0.5 (the fastest tried).
- Turning level keeps up to 0.3 rad/s at 2 mT and 0.5 rad/s at 3-5 mT. At full speed that is a
  turning radius of 13-22 mm.
- Standing up off the floor means levering the whole weight about the tail: 1.02e-5 N m (less
  buoyancy's share). That is 63% of m B / 2, the most a 5 mT field turning across the swimmer
  gives averaged over a turn. A still field does it cleanly:
  - Pointing up, the still field rolls the magnet to the top.
  - Turned slowly back over the swimmer, it levers the nose up with the magnet's full m B.
  - 60 deg of turn stands the swimmer to 39 deg at 5 mT, 33 at 4, 25 at 3 and 10 at 2.
  - Without this stand-up the turning field gets the swimmer off the floor too, but first drags
    it along the floor and the insert: 12.7 s of contact against 0.3 s.

**The course.** A printed insert in the jar: a post 8 mm across on the jar's axis, and two fins
from it to the glass, each with a window 18 mm tall from 9 to 23.5 mm out. The windows are
centred 10 mm below and 10 mm above mid-height, so the whole thing is a spiral staircase.

- The swimmer starts lying on the floor.
- The path for its middle laps the post at 14 mm. It is level within 25 deg of each window and
  climbs 32 deg between them, 98.5 mm in all.
- Two cameras see it at 30 frames a second, with 0.05 mm and 0.5 deg of noise, one frame late.
- The controller is pure pursuit: a target 6 mm ahead on the path, at 6 mm/s. Its aim is tilted
  up against sinking, using the model's drag at the swimmer's seen axis, plus a slow vertical
  correction. The field's rate stays under 0.6 of step-out.
- The field comes from the coils: cfsem on a 2 mm grid, within 7e-5 of cfsem direct.
- Why the windows are 18 mm tall: aiming up against its sinking, the swimmer flies tilted, like
  a helicopter. A 28.8 mm stick crossing a fin sweeps 28.8 mm times the sine of its tilt in
  height (12 mm at 25 deg, 20 mm at 45). With 9 mm windows the first runs jammed against the
  fin, tilted up to 65 deg.

In the table, time touching counts after the body has left the floor. The coils' power is the
mean over the run, and the worst pair's warming assumes no cooling.

| run | finished | time, s | middle's distance from the path, rms / most, mm | touching, s | coils, W | worst coil, C |
|---|---|---|---|---|---|---|
| 2 mT | no (6%) | 41 | 2.5 / 4.0 | 30 | 11 | 0.5 |
| 3 mT | no (28%) | 33 | 2.9 / 6.3 | 22 | 21 | 0.7 |
| 4 mT | yes | 24.8 | 3.3 / 5.4 | 0.9 | 37.5 | 0.97 |
| 5 mT | yes | 20.6 | 2.8 / 4.3 | 0.3 | 62.2 | 1.46 |
| 5 mT, a perfectly uniform field | yes | 20.6 | 2.9 / 4.6 | 0.3 | 62.1 | 1.46 |
| 5 mT, plus a loop on the seen axis | yes | 26.2 | 5.8 / 12.4 | 3.2 | 62.5 | 1.91 |
| 5 mT, the field pointed at the target | yes | 29.8 | 5.4 / 8.8 | 0.9 | 60.0 | 1.82 |
| 5 mT, no stand-up | yes | 32.4 | 2.8 / 5.9 | 12.7 | 63.3 | 2.23 |
| 5 mT, 0.25 ms steps | yes | 20.3 | 2.8 / 4.4 | 0.3 | 62.1 | 1.44 |

- **Robust:** over five camera-noise seeds the 5 mT run finishes in 20.5-20.6 s at 2.74-2.77 mm
  rms. At 4 mT, three seeds finish in 24.7-24.8 s.
- **Barely touches:** at 5 mT it brushes a fin and the glass for 0.3 s, pressing in under
  0.2 mm. The time includes 4.8 s of standing up; it leaves the floor at 5.0 s.
- **The coils' unevenness doesn't matter here:** a perfectly uniform field runs the course
  within 0.1 s.
- **Simple beat clever:** a loop on the seen axis (proportional and integral, to take out the
  droop and swing) made runs worse. The axis settles over seconds and swings sideways when
  pushed, and the loop fought it. Pointing the field straight at the target, with no allowance
  for sinking, sags: twice the path error.
- **Tail first:** in an early version, climbing hard near step-out, the swimmer slipped, tumbled
  and carried on tail first. A cross-magnetized swimmer moves along the field's axis either way
  round, and tail first its heavy nose hangs below like a pendulum. The controller now steers
  the way the swimmer is actually going and keeps the field's rate under 0.6 of step-out.

**Contact in sim-core** (`steer wall`, `steer contact`). The first course runs failed oddly. The
swimmer sank 2-4 mm into the glass and once shot through a fin at 1 m/s. Probes trace it to
sim-core 0.9's contact detection, and it shows at an everyday scale too: a 1 kg body of 5 cm
radius, 2 ms steps, default settings.

- A capsule lying on a box sinks to the box's centre, whichever geom is listed first. In
  `collide_capsule_box` the contact normal's sign is the reverse of sphere-box's, so the contact
  pulls the capsule in.
- A capsule on a cylinder's flat end falls straight through, and so does a cylinder on a box, on
  a cylinder's end or on its side. No contact is found.
- Spheres rest on every shape, and capsules and cylinders rest on a plane.

So the swimmer's contact shape is spheres, 143 of them, no farther apart than their radius. Then
contact behaves: swimming into a wall at 5 mT, the swimmer presses in 0.006-0.07 mm.

One weakness of the linearly implicit drag shows. After hitting the wall, the swimmer slides
along it at 9.6 mm/s with 1 ms steps, 7.2 at 0.25 ms, and 6.4 with the plain drag at 20 us
steps. The drag assumed the swimmer moves freely, so friction comes out too weak. Good runs
barely touch (the 0.25 ms run matches the 1 ms one), so 1 ms is used. Logged in
docs/CORTENFORGE.md.

Not yet:

- walls in the fluid model (near a wall drag rises and a corkscrew drifts sideways: Lauga et
  al., Biophys. J. 90, 400 (2006));
- the cameras' own calibration, and delay beyond one frame (and in a body there is no camera:
  see "Seeing it" below);
- the field's slew (the coils' L/R is 3.2-3.6 ms);
- a print of the insert.

Design lessons for `shape`: a lighter nose, or the magnet nearer the middle, would cut the dart
effect, and a shorter swimmer would turn and fit through smaller windows. Climbing alone was too
narrow a score.

Next: `scale`.

## What's still unsolved (in the field, not just here)

- Seeing it. Tracking something cell-sized deep in a body in real time is hard; ultrasound, MRI
  and X-ray each have trade-offs. Our camera works only because the jar is clear.
- Blood flow. Blood in arteries moves far faster than these swimmers, so most work targets
  slow places: fluid-filled cavities, the eye, tissue, the smallest vessels.
- The body fights back. The immune system clears foreign particles, and every material must be
  safe.
- Getting them out. They have to be retrieved or dissolve safely.
- The size of the machine. Coils that make a strong, even rotating field over a whole torso are
  large and power-hungry (`scale` puts numbers on it).
- Proof. Going from mouse studies to human trials takes years of safety work.

## How we talk about it

- Lead with what it is for, then show how it works, then say what doesn't work yet.
- Plain words, no hype ("could", never "will cure"). Every number from a run, every claim with
  a source.
- At any demo, keep the small neodymium magnets away from kids and pets: swallowing two or more
  can pinch the gut between them.
