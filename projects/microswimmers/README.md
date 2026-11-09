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

`cargo run --release -p microswimmers -- <mode>`, modes added as they are built.

1. `check`: the physics before anything is built on it.
   - Drag: a sphere against Stokes (6 pi mu a v); a slender rod against slender-body theory
     (Batchelor 1970; Cox 1970).
   - The swimmer model: resistive force theory (Gray and Hancock, J. Exp. Biol. 32, 802 (1955);
     Lighthill, SIAM Rev. 18, 161 (1976)) in plain Rust. Its propulsion matrix is symmetric; a
     straight rod spun about its axis doesn't swim; a mirrored helix swims backwards; below
     step-out, speed is linear in the rotation rate.
   - A second model to check it: regularized Stokeslets (Cortez, Fauci and Medovikov, Phys.
     Fluids 17, 031504 (2005)). Resistive force theory drifts for thick helices; map where the
     two agree and use that range.
   - CortenForge's own fluid model. Reading `cortenforge-sim-core` 0.9
     `src/forward/passive.rs` (`mj_ellipsoid_fluid`, `mj_inertia_box_fluid`), the viscous drag
     on a geom is the same in every direction, and the direction-dependent terms are the
     fast-flow (density) ones. So a helix built from sim-core geoms should make no thrust in
     syrup. Confirm by a run; if it holds, log it in docs/CORTENFORGE.md with a Suggestion
     (direction-dependent viscous drag per geom would make resistive force theory native). The
     swimmer's fluid forces then come from our model through `xfrc_applied`.
   - cfsem: reuse coil-designer's Helmholtz checks for the coil field.
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
