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
is the ratio of inertia to stickiness, the Reynolds number. A printed corkscrew a few
millimetres long in glycerol or corn syrup is in the same regime as a microrobot in water (both
far under 1; `swim` computes the two side by side). So:

- a printed helix a few mm long with a tiny magnet magnetized across its axis,
- a jar of glycerol inside three pairs of coils (x, y, z),
- two pairs driven 90 degrees apart make the rotating field; the third tilts it to steer,
- a camera above, to measure speed against the model and to steer through a printed maze.

Something people can watch and understand is easier to trust than a promise.

## Plan

`cargo run --release -p microswimmers -- <mode>`, modes added as they are built (so far:
`check`).

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
2. `swim`: speed against rotation rate for a printable helix in glycerol, up to the step-out
   frequency, where the magnet can't twist it fast enough and it slips and slows. The Reynolds
   numbers of the desk demo and a cell-sized swimmer in water side by side.
3. `coils`: three Helmholtz pairs round the jar (sized with coil-designer): how even the
   rotating field is over the jar, inductance, the drive voltage at the rotation rates `swim`
   needs, power and heating, and the leftover gradient's pull against the twist.
4. `shape`: the helix (radius, pitch, turns, wire thickness) optimized with sim-opt for speed or
   step-out, under printability (cf-design `pipe_spline`, mesh-printability). Writes the STL.
5. `steer`: sim-core with gravity, buoyancy and contact with the jar and a printed maze; a
   controller tilts the rotation plane to follow a path (sim-rl if a hand-tuned one isn't
   enough).
6. `scale`: the honest numbers. Coils sized for a human torso, their power and stored energy,
   and the swimmer's speed against blood flow. This fills the section below with real figures.
7. Hardware: print the swimmer and the coil formers, a three-channel driver, the jar and a
   camera; measure speed against rotation rate and compare with `swim`.

Connects to coil-designer (coils and formers), magnet-tracking (open question: can
magnetometers find the swimmer's magnet if the drive field pauses while they read?) and
colon-probe (the same steer-from-outside idea at a larger size).

## Checks (`check`, 2026-11-09)

23 of 23 pass, about 65 s. In plain words: our model of a swimmer in syrup reproduces exact
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

c = 0.3 has the least worst error at every spacing and is used from here on. (h = 0.785 a is a
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

Next: `swim`.

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
