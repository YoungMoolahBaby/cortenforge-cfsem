# Active radiation shield for the mothership

A magnetic shield for the torus habitat in `C:\dev\mothership` (10,000 people, a 70 m tube whose
centre circle is 824 m from the spin axis, spun at 1 rpm). The shield today is 8.0 t/m^2 of water
for 10 mSv/yr, 18.6 million t, about 95% of the ship. Here: a confined toroidal field between two
windings round the tube (inner at 75 m, outer at 75 + d m, opposite currents), so the field
B = mu0 N I / (2 pi R) lives only in the shell and the habitat is field-free; cosmic rays must
cross it. How much it lets through, what it weighs, and what it needs.

- cfsem: `flux_density_circular_filament_cartesian_par` and
  `vector_potential_circular_filament_cartesian_par` (oriented loops, finite wire radius) for
  the windings' field, stored energy (flux linkage) and forces (I dl x B).
- Ours: the shield geometry (`src/shield.rs`), a Boris tracer in path length with backtracing
  by Liouville's theorem and sphere tracing through the field-free regions (`src/trace.rs`), a
  free-space GCR dose model (`src/dose.rs`), the sizing (`src/size.rs`).
- CortenForge: nothing yet; the magnet's support structure wants a linear-elastic solver it
  doesn't have (logged in `docs/CORTENFORGE.md` Suggestions).

## Modes

`cargo run --release -p active-shield -- <mode>`

- `check` (default): 15 checks against closed forms (about 4 s). Exits 1 on any failure.
- `transmission`: the cosmic-ray transmission against rigidity for the baseline (100 T m, 5 T,
  d = 20 m) and how the cutoff moves with bending power and thickness (about 1 min).
- `size`: conductor, structure, dose through the magnet, water behind it for 10 mSv/yr, and
  cryocooler power for eight designs (about 1 min).
- `failsafe`: water thick enough for a dose limit with the magnet off, and the smallest magnet
  (40 m shell) that brings it to 10 mSv/yr with it on, for limits of 20, 50 and 100 mSv/yr
  (about 1 min). The mothership's design (pass 6) is the 50 mSv/yr row.

## Check (2026-11-06, 15 of 15)

- One loop's centre field is mu0 I / 2a along its normal (3e-12): cfsem's sign convention.
- 2,000 loops (two windings of 1,000) against the continuous winding's mu0 N I / (2 pi R):
  3.4e-5 mid-shell, at a coil and between coils; inside the inner winding the field cancels to
  1.3e-10 of the shell field.
- Stored energy (1/2) sum I Phi, with Phi the line integral of A round a loop and a wire radius of
  the strip's geometric mean distance (0.2231 x the loop spacing), against the closed form
  (mu0 (N I)^2 / 2) (sqrt(R^2 - r1^2) - sqrt(R^2 - r2^2)): 1.3%.
- Magnetic pressure B^2 / 2 mu0 pushes the inner winding in and the outer out: the radial force
  round each loop matches mu0 (N I)^2 r / (2 N sqrt(R^2 - r^2)) to 0.7% and 0.6%.
- The coil pair's pull toward the spin axis, a small difference of large forces, converges to
  (1/N) dE/dR at order 1.00 in the loop count (8.2% / 4.1% / 2.0% at 1,000 / 2,000 / 4,000).
- Boris tracer: a 1 GV proton in 2 T turns on 2 B rho / B to 2e-12, closes after one turn to
  5e-8, keeps its speed to 9e-15. A flat slab (a shell 1,000 km round, far from the axis) stops
  normal incidence below c B d and passes it above (0.97 / 1.03 x 29.9 GV).
- Dose model: proton stopping at 100 MeV 7.232 MeV cm^2/g against NIST PSTAR's 7.289 (0.8%);
  minimum ionisation 2.034 against PDG's 1.981 (2.7%, no density correction); free-space skin
  dose equivalent at solar minimum 129.4 cSv/yr against Ambroglini et al. 2016's 117.7 (9.9%),
  mean quality factor 6.9, iron-nickel the largest group (41.9 cSv/yr).
- Leakage at the habitat's surface, 5 m inside the inner winding: 18 mT with 1,000 loops (4.7 m
  apart there), 0.068 mT with 2,000 (2.4 m), 4.6e-6 mT with 4,000. Living spaces should stay
  under 0.5 mT ("much lower restriction levels, such as 0.5 mT", ICNIRP 2009, for implants), so
  the winding must be ~2,000 turns-equivalent or finer; the Earth's field is 0.03-0.07 mT.

## Transmission (`transmission`, 24 poloidal points x 400 cosine-weighted directions)

The habitat sees 0.979 of the sky with the field off (the ring shades its own inner side: 0.925
at the inner equator); transmissions below are relative to that. Baseline, 100 T m (c B d = 30
GV):

| GV | 1 | 5 | 10 | 15 | 20 | 25 | 30 | 40 | 50 | 70 | 100 | 150 | 200 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| overall | 0.0002 | 0.0012 | 0.0047 | 0.017 | 0.160 | 0.325 | 0.451 | 0.626 | 0.733 | 0.858 | 0.953 | 0.993 | 0.998 |
| inner equator | 0 | 0.0027 | 0.0054 | 0.0081 | 0.089 | 0.260 | 0.384 | 0.576 | 0.670 | 0.808 | 0.922 | 0.981 | 0.984 |

- Half get through at 32.8 GV, 10% at 17.9, 90% at 83.2: the cutoff is soft, and the inner side
  of the ring (stronger field, B ~ 1/R) stops more.
- The half-transmission rigidity scales with bending power, ~1.1 x c B d: 8.3, 16.4, 32.8, 65.5
  GV at 25, 50, 100, 200 T m (d = 20 m).
- Thicker, weaker shells win at the same bending power: 100 T m at d = 10 / 20 / 40 m stores
  1.04e15 / 5.53e14 / 3.10e14 J and lets 90% through at 103 / 83 / 62 GV (sharper).
- A trickle gets through even at 1 GV (2e-4): particles moving along the toroidal field drift
  across it (curvature and gradient drift, ~ r_L / R a metre) over several laps of the ring. A
  pure toroidal field can't hold particles that run along it, the reason tokamaks add a poloidal
  field.

## Size (`size`)

Inputs (sources below): REBCO J_e 1,000 A/mm^2 at 8,870 kg/m^3; carbon-fibre structure 1,200
MPa at 1,550 kg/m^3 with Q = 3; 0.5 t/m^2 of water always inside the magnet, over the habitat's
2.326e6 m^2; heat leak 1.25 W/m^2 onto the inner winding at 43 W/W of cryocooler input; solar
minimum (phi 0.42 GV); secondaries x2 on the primaries-only dose; water's e-fold 2.01 t/m^2.

| T m | d m | B T | 50% GV | energy J | magnet Mt | dose mSv/yr (primaries, x2) | water t/m^2 | hybrid Mt | t/m^2 |
|---|---|---|---|---|---|---|---|---|---|
| 50 | 20 | 2.50 | 16.4 | 1.38e14 | 0.63 | 40, 79 | 4.16 | 10.32 | 4.44 |
| 100 | 20 | 5.00 | 32.8 | 5.53e14 | 2.34 | 14, 28 | 2.03 | 7.07 | 3.04 |
| 100 | 40 | 2.50 | 29.4 | 3.10e14 | 1.31 | 17, 34 | 2.48 | 7.09 | 3.05 |
| 150 | 40 | 3.75 | 44.5 | 6.97e14 | 2.86 | 9, 19 | 1.24 | 5.76 | 2.48 |
| 200 | 40 | 5.00 | 60.1 | 1.24e15 | 5.02 | 6, 12 | 0.50 | 6.18 | 2.66 |
| 200 | 60 | 3.33 | 55.5 | 9.14e14 | 3.70 | 8, 15 | 0.83 | 5.64 | 2.43 |
| 300 | 60 | 5.00 | 83.8 | 2.06e15 | 8.21 | 4, 8 | 0.50 | 9.37 | 4.03 |
| 400 | 80 | 5.00 | 103.4 | 3.01e15 | 11.93 | 3, 7 | 0.50 | 13.09 | 5.63 |

- **A hybrid is ~3x lighter than water alone**: the best, 150-200 T m with d = 40-60 m plus
  ~0.8-1.2 t/m^2 of water, is 5.6-5.8 Mt against 18.6 Mt of water for the same 10 mSv/yr.
  The structure is most of the magnet (the conductor is 0.1-0.26 Mt); it scales with Q: if the
  windings carry the pressure as hoops only, Q is ~9.6 at 75/95 m and 5.3 at 75/125 m (the
  pass-6 researcher's derivation), which would roughly double to triple the magnet.
- **The magnet alone can't reach 10 mSv/yr** at sane sizes once secondaries count: 300-400 T m
  get to 7-8. Published magnets (4-23 T m) reach about 10x below free space (SR2S: "~10 cSv/y"
  at 23 T m); this needs ~100x, so it pairs with water.
- **Cryocoolers: ~131 MW** for the inner winding's heat leak (1.25 W/m^2 over 2.44e6 m^2 at
  ~40 K), as much as the ship uses today; it doesn't change with bending power. Better
  insulation, a warmer intercept shield or warmer HTS operation are the levers.
- **Stored energy 3-12e14 J, 75-290 kt of TNT**: a quench must be caught and the energy spread
  or dumped; the windings would be many independent sections.
- **No braking of the spinning hull**: the field is axisymmetric, so it is static in the hull's
  spinning frame; only the ripple between coils (16.7 Hz at 1,000 loops) induces currents, and
  it is gone at 2,000+. Aluminium rather than steel for the hull near the magnet.
- Solar particle events (protons below a few GV) are stopped outright: no storm shelter needed.

## Fail-safe hybrid (`failsafe`)

The magnet as a bonus layer, not a lifeline: the water alone keeps the dose under a limit if the
magnet or its cryocoolers fail (Globus & Strout's OLTARIS water curve from its 7 t/m^2 row at
2.01 t/m^2 an e-fold); the magnet's leftovers (x2 for secondaries) attenuated by that water meet
10 mSv/yr with it on, never above the magnet-off dose.

| magnet-off limit | water t/m^2 | water Mt | magnet | magnet Mt | shield Mt | on mSv/yr | stored J | kt TNT | cryo MW, 40 / 65 K |
|---|---|---|---|---|---|---|---|---|---|
| 20 mSv/yr | 6.56 | 15.27 | 30 T m, 0.75 T | 0.141 | 15.41 | 7.6 | 2.79e13 | 7 | 131 / 55 |
| **50 mSv/yr** | **4.72** | **10.98** | **50 T m, 1.25 T** | **0.354** | **11.34** | **9.1** | **7.74e13** | **18** | **131 / 55** |
| 100 mSv/yr | 3.33 | 7.74 | 100 T m, 2.50 T | 1.309 | 9.05 | 6.6 | 3.10e14 | 74 | 131 / 55 |

- 50 mSv/yr is ICRP's occupational ceiling for any single year ("no single year exceeding 50
  mSv", ICRP 103 Table 6), so any repair time is survivable; 4.7 t/m^2 of water stops solar
  particle events by itself. The mothership's design: 11.3-11.5 Mt against 18.5 Mt of water
  alone (7 Mt saved), a quarter of the full hybrid's stored energy.
- These magnets are low-field (0.75-2.5 T), so the conductor can run warm, ~65 K (J_e ~300
  A/mm^2 instead of 1,000, judgement from Jc ~ exp(-T/T*)): about 3x the conductor (0.48 instead
  of 0.35 Mt at 50 T m) for 55 MW of cryocoolers instead of 131.
- Magnet leftovers through the 40 m shells: 399, 171, 100, 66, 48, 37, 27, 17, 12, 9 mSv/yr of
  primaries at 10, 20, 30, 40, 50, 60, 75, 100, 125, 150 T m.

## Inputs and sources

- GCR: BPH00 local interstellar spectrum and the force-field approximation (Usoskin et al. 2005,
  eqs. 1-2: https://cc.oulu.fi/~usoskin/personal/2005JA0112501.pdf), phi = 420 MV at solar
  minimum (its Table 3: 1987 420, 1996 440 MV); element groups scaled by PDG 2020 Table 30.1 at
  10.6 GeV/nucleon (https://pdg.lbl.gov/2020/reviews/rpp2020-rev-cosmic-rays.pdf).
- Dose: Bethe (PDG 2024 eq. 34.5, K = 0.307075 MeV mol^-1 cm^2; water I = 79.7 eV, Z/A =
  0.55509), ICRP 60 Q(L); free-space reference Ambroglini, Battiston & Burger 2016, Table 1
  (https://www.frontiersin.org/articles/10.3389/fonc.2016.00097/full): skin 117.7 cSv/yr, BFO
  62.4 at solar minimum. Primaries only; the x2 for secondaries is a judgement.
- Conductor: "Record JE values for commercial wires of over 1000 A/mm2 at 20 K, 20 T" and "JE
  is in the 500-1400 A/mm2 range" (Molodyk et al. 2021, Sci. Rep. 11:2084); tape density from
  its layers (40 um Hastelloy, 10 um Cu, REBCO, Ag).
- Structure: Toray T1100G composite 3,460 MPa tension (data sheet), 1,200 MPa allowable
  (judgement, safety factor ~2-3); virial theorem M >= rho E / sigma ("if some structure is
  under compression then additional structure must be in tension to compensate", Hassenzahl
  1988); Q = 3 with radial ties across the shell (judgement).
- Cooling: JWST's instruments "cooled to ~40 K ... passively by a tennis court-sized sun shield";
  cryocoolers at 10-20% of Carnot at 40 K (ter Brake & Wiegerinck 2002); heat leak through
  multilayer insulation at an effective emissivity ~0.003 from the 293 K habitat (judgement).
- People: ICNIRP 2009 static-field guidelines (400 mT public, 0.5 mT for implants).
- Water's e-fold of dose deep in water (2.01 t/m^2 from Globus & Strout 2017 Table 2, 6-7 t/m^2,
  where the spectrum is hard like the magnet's leftovers).

## Open

- Transport (secondaries in the water, hull and body) for the transmitted GeV spectrum, in
  place of the x2: an OLTARIS or Geant4 run.
- The structure itself: stresses and Q for real windings and ties (needs a linear-elastic solver).
- Field-aligned drift leakage over longer paths (here capped at 60 km).
- Quench protection, the cryogenic heat budget, and the winding's ripple at the hull.
