# Wireless charging (main 4)

Inductive power for the maglev cart first: a loop along the track, a pickup coil on the cart, so
the cart needs no battery swaps and has no sliding contacts (no brush dust, which suits the
hospital). The same models then serve a charging pad and receiver (and could power the capsule
in colon-probe).

- cfsem: self- and mutual inductance of the track loop and the pickup (straight filaments),
  coupling over the cart's sideways and vertical motion, stray field, induced voltage in the lift
  coils, the steel rail's effect bounded by images.
- CortenForge: a printable pickup former and track tray (cf-design), and the pickup on the cart in
  the maglev simulation (sim-core).

Run: `cargo run --release -p wireless-charging -- <mode>`.

- `check`: cfsem's straight-filament inductance and field against textbook results (about 1 s).
- `link`: the track-to-cart link: a check against the 2D two-wire flux, spacing, height and turns
  sweeps, the chosen design, and its tolerances (about 8 s).
- `effects`: the steel rail at 85 kHz (a 2D conducting magnetic half-space, checked against the
  images), induced voltage in the lift and stator coils, stray field against ICNIRP (about 23 s).

## Checks (`check`, 2026-11-02)

9 of 9 pass.

| Check | Result |
|---|---|
| Rectangle 100 x 50 mm, 1 mm wire, self-inductance | cfsem 277.90 nH; Grover (partial inductances, DC) 263.10 nH: +5.6%. Less mu0/8pi per metre: -7.6e-4 |
| Ring 50 mm radius, 1 mm wire (1024-gon) | cfsem 325.76 nH; Wien 310.06 nH: +5.1%. Corrected: -1.1e-5 |
| Coaxial rings 50 and 40 mm, 20 mm apart, as polygons vs cfsem's circular filament | 16-gons -2.2e-2, 64 -1.4e-3, 256 -8.8e-5, 1024 -5.5e-6 (second order) |
| Two offset squares (3 and 2 turns) vs own Gauss-Legendre Neumann integral | +2.1e-10, reversed +1.4e-9 |
| 20 m track loop field vs finite-wire Biot-Savart | to 1e-15 (the 2D two-wire limit is off by up to 4e-5 from the 10 m ends) |
| Images across a plane: conductor (opposite current), iron (same current) | normal B 2.3e-16, tangential B 1.3e-16 of the peak |

- cfsem's linear-filament self-inductance is about 5% high for single turns of 1 mm wire: it
  evaluates A on the wire's axis, which counts the internal inductance as mu0/4pi per metre
  instead of mu0/8pi (logged in docs/CFSEM.md). The excess is exactly mu0/8pi x length, so
  `loops::self_inductance` subtracts it. Mutual inductances are unaffected.
- My first image helper reversed the iron image's loop order (the wrong current sense); check 6
  caught it (tangential B 1.3 x the peak).

## The track-to-cart link (`link`, 2026-11-02)

A one-turn Litz track loop in a tray under the cart (two conductors 60 mm apart, from x = -0.15
to 1.35 m, so past both ends of the 0 to 1.2 m trip) and a flat 10-turn Litz spiral pickup, 130 x
70 mm, under the cart's hull, 8 mm over the track conductors at the hover gap. 85 kHz (the
SAE J2954 band). The track is driven at constant current through a series capacitor; the pickup
is parallel-tuned, which raises its voltage by Q2 to suit a DC bus. Assumed (to replace with
datasheets): Litz AC/DC resistance 1.2, a resistive load standing for rectifier and charger.

- One pickup turn over the middle of a 20 m track matches the 2D two-wire flux times its length
  to 4.5e-6 (at 40, 60 and 80 mm spacing); on the real 1.5 m track the ends add about 1e-3.
- Pieces of 5 mm vs 2.5 mm on the pickup change L2 by 2e-5 and M not at all (to 4 digits).
- Track loss over load power is 1 / (k^2 Q1 Q2), with k^2 Q1 = w M^2 / (L2 R1), which no choice
  of turns changes. Spacing sweep (12-turn pickup, 8 mm under it): k^2 Q1 0.30 at 30 mm, 0.97 at
  60 mm, 0.54 at 90 mm, so 60 mm (the pickup's mean turn width); best coil-to-coil efficiency
  70.3%, 80.9%, 76.0%.
- Height: 5 mm under the pickup 83.4%, 8 mm 80.9%, 12 mm 77.7%, 16 mm 74.5%. Chosen 8 mm: room
  for a tray cover, and the cart only ever rises from there (1.4 mm on landing).

Turns set the voltage, at loaded Q2 = 8 and 5 W into the load:

| Pickup turns | L2 | M | R_L | Load voltage | Track current | Coil to coil |
|---|---|---|---|---|---|---|
| 6 | 8.38 uH | 0.594 uH | 35.8 ohm | 13.4 V rms | 6.15 A | 73.3% |
| 8 | 13.16 uH | 0.769 uH | 56.2 ohm | 16.8 V rms | 5.81 A | 76.1% |
| 10 (chosen) | 18.37 uH | 0.923 uH | 78.5 ohm | 19.8 V rms | 5.63 A | 77.8% |
| 12 | 23.75 uH | 1.055 uH | 101.5 ohm | 22.5 V rms | 5.54 A | 78.8% |
| 14 | 29.13 uH | 1.165 uH | 124.5 ohm | 24.9 V rms | 5.52 A | 79.4% |

Chosen 10 turns (3.57 m of Litz): 19.8 V rms rectifies to about 17 V for the 12 V battery's
charger. C2 190.9 nF (19.8 V rms, 2.04 A rms in the coil). Track: 5.63 A rms, series C1 1.37 uF,
7.7 V rms across its 2.57 uH. Losses: track 0.82 W, pickup 0.61 W.

With the capacitor, load and track current fixed:

| Case | M | L2 | Load power |
|---|---|---|---|
| middle of the track | | | 5.00 W |
| at either end of the trip (x 0 or 1.2 m) | +1.11% | | 5.11 W |
| 2 mm sideways (at the lips) | -1.05% | | 4.90 W |
| landed on the skids, gap 4.53 mm | -7.26% | | 4.30 W |
| payload on, mean gap 5.55 mm | -2.20% | | 4.78 W |
| on the station shelf, gap 6.5 mm | +2.98% | | 5.30 W |
| rail as ideal iron | +10.64% | +5.30% | 4.98 W |
| rail as a perfect conductor | -10.64% | -5.30% | 3.76 W |

- The cart's own motion barely matters: along the whole trip and up to the lips, within +-2%.
- The steel rail matters most, 38 mm above the pickup. At 85 kHz steel is both magnetic and a
  good conductor with a skin depth of tens of um, so it sits between the two bounds, and it will
  also heat (eddy currents): see Side effects below.
  Sizing for the conductor bound means about 15% more track current; the steel model gives 6.61 A.

## Side effects (`effects`, 2026-11-02)

### The steel rail at 85 kHz

cfsem's images give the bounds; real steel is magnetic and conducting with a skin depth of
20-70 um. Over a conducting magnetic half-space each spatial frequency k across the track is
reflected with weight R(k) = (mu_r k - gamma) / (mu_r k + gamma), gamma = sqrt(k^2 + j w mu0 mu_r
sigma) (R = 1 is ideal iron, -1 a perfect conductor). For the long track this is a 2D integral;
its ratio to the iron case scales cfsem's 3D iron image terms, and its imaginary part is eddy loss.
Steel: sigma 6 MS/m (assumed), mu_r 100 to 1000 (uncertain; the rail carries the magnets' bias
flux).

- Checks: the spectral integral with R = 1 matches the closed-form image to 5e-8 (track to
  pickup, track and pickup to themselves), and R = -1 its negative to 5e-8; the 2D image times
  each turn's length gives 0.0975 uH for the track-to-pickup term against cfsem's 3D 0.0982 uH.
- The rail acts mostly as a conductor:

| mu_r | Skin depth | Image weight (M) | M | L2 | Track for 5 W | Rail heating | Track to load |
|---|---|---|---|---|---|---|---|
| 100 | 70 um | -0.865 - 0.114j | -9.20% | -4.53% | 6.58 A | 0.56 W | 68.7% |
| 300 | 41 um | -0.773 - 0.173j | -8.22% | -4.01% | 6.61 A | 0.84 W | 66.0% |
| 1000 | 22 um | -0.612 - 0.252j | -6.51% | -3.12% | 6.61 A | 1.20 W | 63.0% |

- Sized for the worst case: 6.61 A rms in the track (stage 2's 5.63 A assumed no rail). The rail
  costs 6-15 points of efficiency and warms by 0.6-1.2 W over the track's length.
- The 2D ratio is also applied to the pickup's own image (L2), which is less exact: the pickup is
  130 mm long, not long against its 76 mm distance to its image.
- Not modelled: a ferrite backing on the pickup, which real pickups use to shield and to raise
  coupling. cfsem has no magnetic materials (only infinite-plane images, which would overstate a
  130 x 70 mm sheet), so it needs another solver.

### Induced 85 kHz voltage

At 6.61 A in the track and 2.39 A rms in the pickup (-96 deg from the track current), cart at
mid-track; the three values are the rail free / ideal iron / perfect conductor:

| Coil | Voltage (rms) | Notes |
|---|---|---|
| each lift coil (80 turns) | 702 / 1125 / 279 mV | 10.8 / 17.3 / 4.3 mA through its 0.122 mH |
| guidance coil under the cart (20 turns) | 964 / 1854 / 75 mV | steel model: 277-604 mV |
| guidance coil 0.3 m from the cart | 185 / 357 / 13 mV | steel model: 53-116 mV |
| motor coil under the cart (20 turns) | 330 / 648 / 12 mV | steel model: 72-173 mV |

- The four lift coils see the same magnitude (symmetry). 4-17 mA at 85 kHz is the same size as the
  hovering control currents (about 13 mA), far above the 1 kHz current loop but sampled at 2 kHz:
  the current and gap measurements need anti-alias filtering ahead of the sampling, or the
  ripple folds into the control band. The force it makes averages out at 85 kHz.
- Stator coils lie on the rail's face, where the rail sets the normal field (iron doubles it, a
  conductor cancels it); the steel model gives 0.22-0.63 of the free value. Under a volt per coil:
  the drivers must tolerate it; the circulating current and its loss were not computed.
- Gap and sideways sensors (eddy-current class) work at their own carrier frequency and need a
  check against an 85 kHz field; not modelled.

### Stray field

ICNIRP 2010 reference levels (Health Physics 99(6), Tables 3 and 4), 3 kHz-10 MHz, unperturbed
rms: 27 uT (21 A/m) for the general public, 100 uT (80 A/m) occupational. At 6.61 A, worst of
the rail bounds, distances from the track conductors:

| Where | 5 cm | 10 cm | 20 cm | 50 cm | Over 100 uT to | Over 27 uT to |
|---|---|---|---|---|---|---|
| below the track, cart over it | 85 uT | 25.0 uT | 5.4 uT | 0.69 uT | 4 cm | 9 cm |
| beside the track, cart over it | 39 uT | 12.0 uT | 3.3 uT | 0.48 uT | 2 cm | 6 cm |
| below the track, 0.4 m from the cart | 27 uT | 9.7 uT | 3.1 uT | 0.58 uT | 0 cm | 5 cm |
| beside the track, 0.4 m from the cart | 16 uT | 6.3 uT | 2.4 uT | 0.40 uT | 1 cm | 3 cm |

- The go-and-return pair cancels fast (falls as 1/r^2): over the public level only within 9 cm of
  the conductors, so a tray enclosure 10 cm deep keeps the public level outside it, unaveraged.
- Medical implants (pacemakers and similar) have their own limits, not checked here; a hospital
  would need that before installing.

## Plan

1. Checks (done).
2. Track loop and pickup for the cart (done): the 85 kHz link, coupling over the trip and the
   cart's motion, the steel rail bounded by images.
3. Side effects (done): the steel rail (a conductor, mostly; 6.61 A for 5 W, 0.6-1.2 W of rail
   heating), induced voltage in the lift and stator coils (anti-alias filtering needed), stray
   field (public level within 9 cm).
4. Printable pickup former and track tray; the pickup on the cart in the maglev simulation.
