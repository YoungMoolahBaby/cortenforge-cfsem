# Maglev delivery cart, hybrid EMS (main 2)

A contact-free hospital delivery cart: the kind of small track vehicle that carries lab samples,
blood and medicines between wards, but levitated, so nothing rubs, wears, sheds or needs
grease. Built as a printable tabletop demonstrator.

- cfsem: magnet and coil fields and forces against the steel rail (iron by images), coil
  inductance for the current dynamics, stray-field maps, and later the linear-motor thrust.
- CortenForge: cart and track simulation (sim-core), the gap controller, printable track and
  cart parts (cf-design).

## Why hybrid EMS

Electromagnetic suspension (EMS) pulls the cart up towards a steel rail and holds the gap by
feedback. With air-core coils alone the pull is too weak at this scale (`check`: a 408-turn coil
makes 0.049 N at 1 A and 5 mm, so a 300 g cart would burn ~170 W). Real EMS uses iron-core
magnets, which cfsem cannot model (it has no magnetic materials). The hybrid puts permanent
magnets under the rail to carry the weight and small coils around them to stabilize: a 20 x 5 mm
N42 disc pulls 8.9 N at 5 mm, and a 72-turn coil adds 0.30 N per amp. A well-tuned controller
drives the mean coil current to zero, so the magnets carry the load for free ("zero-power"
control). Magnets are modelled as their equivalent surface currents, which is close to exact for
NdFeB (mu_r about 1.05).

## Reference: hospital track systems (electric track vehicles)

Today's systems run on wheels along wall and ceiling tracks, powered through conductor rails,
climbing on a rack; vendors include Telelift and Swisslog. Published figures vary by source:
speeds up to 1 m/s, payloads around 10-15 kg for the small cars and 40-50 kg for the large ones,
horizontal and vertical travel in one vehicle; cargo is lab samples, blood, pharmaceuticals,
sterile supplies and records, mostly under 10 kg. Sources: [ETV systems (Wikipedia)](https://en.wikipedia.org/wiki/Electric_track_vehicle_system),
[Telelift (Wikipedia)](https://en.wikipedia.org/wiki/Telelift), [Transtech ETV](https://transtech-sys.com/product/etv/?lang=en),
[Practical Patient Care](https://practical-patient-care.com/?p=4803).

## Requirements for an ultra-clean cart

For a strict hospital (think Mayo Clinic), at full scale, and scaled down for the demonstrator:

- No contact in normal running: no wheel or rail wear particles, no lubricant, no brush dust.
- Smooth, sealed, washable surfaces (hospital disinfectants); no fans to stir air.
- Quiet, low vibration (samples, especially blood).
- Low power at rest and cruise (the magnets carry the load).
- Stray field: below 0.5 mT at any surface a person can reach, the "5 gauss" line used for
  pacemaker safety around MRI scanners. To be mapped with cfsem.
- Safe on power loss: the magnets keep pulling, so the cart must land on skids, not fall.
- Payload and speed: full scale ~10 kg at up to 1 m/s (small ETV class); the demonstrator a few
  hundred grams at walking-pace scale.
- Known gaps for later: vertical runs (ETVs climb on a rack), track switches, the linear motor.

## Modes

`cargo run --release -p maglev-train -- <mode>`

- `check` (default): force and inductance model against known results, 19 checks, plus design
  numbers for air-core vs hybrid. Exits 1 on any failure.
- `stray`: magnet sizing for the demonstrator, whole-cart lift, the 0.5 mT reach for all-same
  and alternating polarity, rail flux, and field maps in `out/stray_*.svg`.
- `corner`: one corner in CortenForge's simulator: open-loop instability vs prediction, gap
  control through the coil, zero-power control, a sensor-noise sweep; plots in
  `out/corner_*.svg` (about 4 s).
- `cart`: the whole cart: lift-off, an off-centre payload, a 1.2 m trip to a station driven by the
  linear motor, against the zero-power attitude predicted from the force table; plot in
  `out/cart.svg` (about 35 s, most of it building the motor tables).
- `motor`: the linear motor: cfsem's circular-to-linear mutual against the true midpoint,
  per-phase thrust, thrust by Lorentz force vs coenergy, the commutated thrust constant and
  ripple, copper loss (about 22 s).
- `parts`: the landing-skid window from the force table, then the printable parts (cf-design),
  checked and written to `out/parts/*.stl` (mm), and a mass budget (about 35 s).
- `guide`: the sideways-guidance coil layer: per-phase sideways force by Lorentz vs coenergy, the
  front and rear blocks' force and yaw per amp, side effects, and the motor's and guidance's
  behaviour 1 mm off-centre (about 16 s).

## The model (`src/em.rs`)

Air-core coils (one cfsem filament per turn) and disc magnets (equivalent current loops) below a
steel rail whose flat underside is z = 0. An ideal soft-iron rail (wide, thick, unsaturated)
acts as a mirror: the field below it is the sources plus their images with the same current
sign. Force on the cart = Lorentz force from the images (F_z = -2 pi r I B_r per loop);
inductance near the rail = L0 + M(coil, image). Limits: no sideways restoring force from a flat
plane (side guidance needed), no eddy currents in the steel yet, rail must stay unsaturated.

## What `check` found (2026-11-02)

- Force between coaxial loops three independent ways (Lorentz with cfsem's B_r, energy with
  cfsem's flux, the derivative of Maxwell's closed form): agree to 1e-8; the close-loop limit
  matches parallel wires (mu0 I^2 a / d) to 3e-6.
- The image across an iron surface: B_r on the plane is exactly zero, B_z is not.
- A 408-turn AWG 24 pancake coil (r 6-19.7 mm, 9.7 mm tall, 2.77 ohm): lift by Lorentz and by
  coenergy (-1/2 I^2 dM/dg) agree to 1e-8 at 2, 5, 10 mm (0.098, 0.049, 0.018 N at 1 A);
  free-space inductance 2.92 mH, 1.6% from Wheeler's multilayer formula.
- Design numbers: air-core, a 300 g cart on 4 coils needs 3.88 A each, 167 W. Hybrid, a 20 x 5 mm
  N42 disc (Br 1.3 T) pulls 25.6 / 8.9 / 3.7 / 1.4 N at 2 / 5 / 8 / 12 mm (the magnet model
  converges: 8.886 / 8.891 / 8.892 N with 10 / 40 / 160 loops); the rail surface sees at most
  0.50 T at 2 mm (steel saturates around 1.5 T). A 72-turn AWG 24 coil around the disc (0.49 ohm)
  adds 0.30 N per amp at 5 mm.

- A disc magnet plus its image as current loops vs two of cfsem's point dipoles: 1.5e-5 to 1.1e-4
  apart at 0.14 to 1 m, inside the O((size/distance)^2) the dipole picture neglects.

## Sizing and stray field (`stray`, 2026-11-02)

Demonstrator: 0.5 kg (cart 300 g, payload 200 g) on four N42 discs at the corners of
100 x 50 mm, 1.226 N a corner.

| Disc (d x h) | Gap carrying a corner | Dipole moment |
|---|---|---|
| 10 x 3 mm | 3.96 mm | 0.244 A m^2 |
| 12 x 4 mm | 5.96 mm | 0.468 A m^2 |
| 15 x 5 mm | 8.94 mm | 0.914 A m^2 |
| 20 x 5 mm | 12.64 mm | 1.625 A m^2 |

Chosen: 12 x 4 mm at 5.96 mm (the smallest stock disc that floats a corner at >= 5 mm).

| Polarity | Cart lift (each magnet vs every image) | vs 4 x single | 0.5 mT reach: down from a magnet / sideways from its edge / 45 deg / down from the cart centre |
|---|---|---|---|
| all same | 4.875 N | -0.61% | 74 / 57 / 71 / 21 mm |
| alternating | 4.934 N | +0.58% | 59 / 46 / 57 / 0 mm |

- Soft iron pulls either pole, so alternating polarity costs no lift (slightly more, here) and
  trims the 0.5 mT reach by about 20% at this size; the cart centre falls below 0.5 mT outright.
  At this size the near field of each magnet dominates; the far-field cancellation matters more
  for a full-size cart with bigger magnets.
- For the demonstrator: reachable surfaces at least ~6 cm from the magnets (alternating), i.e.
  a cart body about that deep below the magnet plane, or a shield.
- A steel keeper or shield under the magnets is the usual way to cut the downward field, but
  cfsem has no magnetic materials, so it can't be checked here (finite iron); noted in
  `docs/CFSEM.md`.
- Rail flux: 7.9e-5 Wb of each sign enters the rail face, 4.0e-5 Wb a magnet; a 60 mm wide rail
  below 1.2 T needs >= 0.55 mm of steel, so any practical flat bar carries it (the image model
  wants it thick anyway).
- Maps: `out/stray_same.svg`, `out/stray_alternating.svg` (x-z plane through a magnet row,
  contours 5 / 1 / 0.5 mT).

## One corner in sim-core (`corner`, 2026-11-02)

The corner: 0.125 kg on the 12 x 4 mm N42 disc, with an 80-turn AWG 28 control coil wound round it
(r 6.5-9.5 mm, 0.859 ohm). Its pull F(g, I) = a(g) + b(g) I + c(g) I^2 is tabulated from cfsem on
2-12 mm (magnet and coil against their images), with the coil inductance L(g) and the magnet flux
Psi(g) it links, which gives the coil equation V = R I + L dI/dt + (I L' + Psi') dg/dt. The body
is a sim-core slide joint (`qfrc_applied` each step), the coil current integrated alongside.

- Linearized at the zero-power gap 5.964 mm: negative stiffness k = 483.3 N/m, open-loop pole
  p = sqrt(k/m) = 62.2 /s (an offset doubles every 11.1 ms); force per amp b = 0.0706 N/A;
  L = 0.122 mH, L/R = 0.142 ms. Energy consistency: b = -dPsi/dg to 1.2e-5, c = -L'/2 to 1.0e-5.
- Open loop, a 10 um offset released from rest (grows as cosh(p t)): 62.04 / 62.02 /s at 50 /
  25 us steps against 62.18 predicted; with the coil shorted, 62.02 / 62.00 against the root
  of L s^3 + R s^2 + (b^2/m - k L/m) s - k R/m = 0, 62.16 (the shorted coil barely damps it:
  b is small).
- Control: a gap PD at 2 kHz (wn = 3p, damping 0.7) commands current (limited to +-10 A); a PI
  current loop at 1 kHz with anti-windup drives the coil from +-12 V; lift-off ramps the setpoint
  over 0.3 s. Zero-power: the setpoint drifts against the mean current (tau 0.5 s) until the
  magnet alone carries the load. The first version, without the current limit, anti-windup and
  ramp, lost the corner even without noise.
- With a 1 um gap sensor and a 100 Hz derivative filter, 4 s runs with +25 g at 1.5 s:

  | Mode | Before payload | After +25 g |
  |---|---|---|
  | gap hold | 5.964 mm, 0.4 um rms, 75 mW | 6.018 mm, 3.71 A, 11.9 W |
  | zero-power | 5.964 mm, 0.4 um rms, 74 mW | 5.512 mm (table: 5.509), +0.02 A, 80 mW |

  Zero-power carries the payload on about 150 times less power.
- The sensor sets the power: coil power from sensor noise alone (zero-power, no payload):

  | Sensor noise | D filter | Gap rms | Current rms | Coil power |
  |---|---|---|---|---|
  | 1 um | 100 Hz | 0.43 um | 0.29 A | 74 mW |
  | 1 um | 300 Hz | 0.39 um | 0.61 A | 323 mW |
  | 2 um | 100 Hz | 0.87 um | 0.59 A | 298 mW |
  | 2 um | 300 Hz | 0.81 um | 1.23 A | 1291 mW |
  | 5 um | 100 Hz | 2.21 um | 1.47 A | 1858 mW |
  | 5 um | 300 Hz | 2.50 um | 3.05 A | 7987 mW |

  Power grows as the square of the noise, because the small force per amp needs high gains.
  Requirement: a gap sensor with about 1 um noise at 2 kHz. For a given winding volume the
  wire choice does not help (b^2/R is fixed); a bigger coil or a model-based estimator would.

## The cart in sim-core (`cart`, 2026-11-02)

One body (140 x 80 x 30 mm, 0.5 kg with its nominal payload; I_xx 3.04e-4, I_yy 8.54e-4 kg m^2)
on six joints: slide along the track, sideways, heave, yaw, pitch, roll. (Before the guidance
layer below existed it had four, with sideways and yaw assumed guided.)
Sites mark the four magnet faces, 15 mm above the centre of mass, at +-50 / +-25 mm. Inside
sim-core's `cb_passive` each corner's gap and gap rate come from its site and `mj_jac_site`, and
the cfsem force table's pull F(g, I) is added to the joint forces as J^T F; the linear motor's
wrench (from `motor`'s tables at the cart's position and mean gap, at the centre of mass) and
the payload's weight go in the same way (`mj_jac`). Coil currents and the controller run each step outside, sharing state
with the callback.

Control is modal. A rigid cart has three vertical motions and four magnets, so four independent
zero-power loops would fight in the fourth pattern (warp), which the body cannot follow. Heave,
pitch and roll each get a gap PD (wn = 3 p of that mode, from its equivalent per-corner mass:
M/4, I_yy/(4 lx^2), I_xx/(4 ly^2); open-loop poles 62.2, 75.2, 63.0 /s) and a zero-power loop;
corner currents are built from those three, so warp current is zero by construction. Same coils,
sensor (1 um, 2 kHz, 100 Hz derivative filter), current loops and limits as `corner`.

Timeline: lift-off ramp to 0.3 s; +100 g set down at (+30, +15) mm at 1 s; depart at 2.5 s for a
station 1.2 m away (0.3 m/s, 0.5 m/s^2): the trip controller asks for a force (mass feedforward
plus PD on a trapezoidal profile) and the linear motor gets that force over its thrust
constant as current. (Stage 3 first ran with a stand-in force at magnet height; stage 4 put the
motor in, with the same results to within 1 um.)

| Phase | Corner gaps (RR / FR / FL / RL) | Corner currents | Coil power |
|---|---|---|---|
| hovering, 0.7-1.0 s | 5.964 / 5.964 / 5.964 / 5.964 mm | within +-0.013 A | 205 mW |
| trip, 2.5-8.0 s (with the 5 s shove) | largest excursion 70.5 um | | 369 mW mean |
| at the station, 8.5-9.0 s | 6.047 / 5.551 / 5.034 / 5.530 mm | within +-0.009 A | 201 mW |

- Stops at x = 1.2000 m. The motor peaks at 0.805 A of thrust current; the stator spends 1.55 J
  on the whole trip (peak 1.56 W). With the payload on, the cart settles tilted (down towards the load) in
  the attitude where the magnets alone carry everything: the currents return to zero.
- That attitude, solved independently from the force table (total force and two moments, three
  mode gaps), matches the simulation to 0.5 um without the payload and 0.7 um with it. The first
  version of that solver missed by 17 um: it left out that the centre of mass sits 15 mm below
  the magnet faces, so tilting moves the weight's line of action (about 6% of the roll moment).
- Setting 100 g down in one instant drives the nearest corner's current to its 10 A limit for a
  moment; a real load set down over a fraction of a second would be gentler.
- Coil power is the sensor noise through the gains (about 50 mW a corner); the trip mean also
  carries the shove's roll transient (sideways force acts 15 mm above the centre of mass). Numbers
  from the six-motion run; the four-motion run before guidance gave 195 / 241 / 225 mW.

## The linear motor (`motor`, 2026-11-02)

The lift magnets alternate polarity along each row, which is the field a linear synchronous
motor needs: the Transrapid trick of letting the levitation magnets double as the motor's field.
A flat three-phase stator lies on the rail's underside in the gap (0.75 mm below the steel; a
printed-circuit coil board, say): full-pitch coils 24 mm across, 20 turns of AWG 24, one coil
track under each magnet row with phases reversed between rows, as the rows' polarities are.
Coils are cfsem linear filaments mirrored in the steel; thrust is the Lorentz force on the
magnets' equivalent loops.

- First design: coil pitch = magnet spacing (100 mm). Each phase's thrust came in narrow spikes
  (195 mN/A when a 12 mm magnet sat right under a coil side, under 10 mN/A elsewhere), so no
  commutation could make it smooth (354% ripple).
- Second design: pole pitch 100/7 = 14.3 mm, close to the magnet size, and an odd fraction of
  the magnet spacing so the two magnets in a row sit pi apart electrically, as their opposite
  polarities need. Per-phase thrust is then sinusoidal (205.9 mN/A peak a phase).
- Thrust, two ways for phase A at 1 A: the Lorentz force on the magnet loops, and the coenergy
  route I dPsi/dx with Psi from cfsem's `mutual_inductance_circular_to_linear` on coil sides
  split into short pieces. With cfsem's function they differ by 9% at 1 mm pieces and 2% at
  0.25 mm; with the segment midpoint taken correctly, by 4.0e-4 and 2.6e-5. That function places
  the midpoint wrongly (40% low on a tangential segment, wrong sign across the -x axis): logged
  as a bug in `docs/CFSEM.md`.
- Commutated (sinusoidal in pi x / tau): thrust constant 0.310 N/A, ripple 0.10% peak to
  peak, extra lift within +-0.16 mN/A, a steady pitch moment of 4.76 mN m/A (the thrust acts
  15 mm above the centre of mass). Copper loss with block switching (only the coils under the
  cart energized): 2.40 W/A^2, 0.129 N/W.

## Sideways guidance (`guide`, `cart`, 2026-11-02)

A flat steel rail pulls straight up only, so with nothing else the cart has no sideways or yaw
restoring force at all: it would drift to a side lip and rub. Real EMS gets lateral force from
iron magnets at the rail's edges, which cfsem can't model; here a second coil layer on the stator
board pushes sideways on the lift magnets themselves (a small x-y planar motor, as in industrial
maglev conveyors):

- Guidance coils 14.3 x 50 mm, 20 turns (0.22 ohm), 1.25 mm under the steel, at the motor's pitch
  and three phases. Their long sides run along the track right over the two magnet rows; the rows
  have opposite polarity and the sides carry opposite currents, so both push the same way. Their
  short sides cross both rows, so their pushes along the track cancel. A plain wire along the
  track would not do: the magnets alternate front to back too, so it gives pure yaw.
- The coils under the cart's front 80 mm and rear 80 mm form two independently driven blocks:
  together they push sideways, opposed they turn the cart.
- Lorentz force vs coenergy (true midpoint) for one phase: 1.2e-3 apart at 1 mm pieces, 7e-5 at
  0.25 mm (second order).
- Per amp of a block: sideways 0.059-0.060 N (2% ripple), yaw +-2.9 to 3.0 mN m; thrust and lift
  about 1e-16 (zero by symmetry); copper 0.93 W/A^2 a block.
- 1 mm off-centre: the motor makes no sideways force to speak of (<= 1e-8 N/A), so the cart is
  neutral sideways; motor thrust falls 0.2%, guidance force 1.7%.

In the cart: sideways sensors at x = +-60 mm, a PD loop on sideways position and yaw (8 Hz,
damping 0.8, on the cart's mass and yaw inertia) giving a force and a moment, solved for the two
block currents from the tables at the cart's position (limit 8 A). Test: start 1 mm off and
0.5 deg turned; a 0.3 N shove at the front edge for 50 ms at 5 s (mid-trip).

| Setting | After the shove (worst body corner) | Contact (lips at 2 mm) | Guidance while hovering |
|---|---|---|---|
| 3 A limit, 4 Hz, 5 um sensors | 2.731 mm | yes | |
| 8 A limit, 8 Hz, 5 um sensors | 0.643 mm | no | 0.707 A rms, 465 mW |
| 8 A limit, 8 Hz, 1 um sensors (chosen) | 0.642 mm | no | 0.141 A rms, 19 mW |

- At 3 A both blocks make at most about 18 mN m of yaw, less than the shove's 21 mN m, so no
  tuning could have held it; 8 A (the lift coils' driver class) does.
- The start offset settles to 0.1 um and 0.0001 deg rms; the largest excursion of the run is the
  start itself, 1.611 mm (1 mm plus 70 mm x sin 0.5 deg).
- Sideways force acts at the magnets, 15 mm above the centre of mass, so it also rolls the cart a
  little: with 5 um sideways sensors the hovering lift coils rose to 250 mW; with 1 um sensors they
  are at 205 mW. Requirement: sideways sensors of about 1 um, the gap sensors' class.
- The rest of the trip is unchanged: stop at 1.2000 m, zero-power attitude within 0.7 um.

## Parts and the parts list (`parts`, 2026-11-02)

### Landing skids

On power loss the magnets keep pulling, so the cart rises until skids on its deck touch the
stator board. They must hold it far enough off that the coils can push it back off (the magnets
pull harder the closer they are), yet not touch in normal running. Push-off current for an
unloaded corner (the hardest case) from the force table:

| Landed gap | 3.5 mm | 4.0 mm | 4.5 mm | 5.0 mm | 5.5 mm |
|---|---|---|---|---|---|
| Push-off current (limit 10 A) | 12.83 A | 10.73 A | 8.43 A | 5.87 A | 3.00 A |

Push-off within 5 A needs a landed gap of 5.16 mm or more, but with 100 g off-centre a corner runs
at 5.03 mm (`cart`): no room at that margin. Chosen: land 0.5 mm under the lowest running gap,
4.53 mm (skids 3.03 mm tall over the deck; the board is 1.5 mm), which needs 8.27 A of the
10 A to push off. Ways to widen the window: a centred or lighter payload (the off-centre load
is what drives one corner down), or a bigger control coil.

### Printed parts

Modelled in cf-design (mm), meshed by marching cubes, welded, checked with the FDM default as
printed (z up), written to `out/parts/`:

| Part | Printed | Volume (solid PLA) | Check |
|---|---|---|---|
| `cart_deck` | flat, pockets up | 61.5 cm^3 (76 g) | printable |
| `cart_hull` | open side up | 51.0 cm^3 (63 g) | printable |
| `rail_hanger` (100 mm) | on end | 79.5 cm^3 (99 g) | 2 flags, see below |
| `station_hanger` (100 mm) | on end | 88.5 cm^3 (110 g) | 2 flags, see below |

- Deck: 140 x 80 x 6 mm, four 19.5 x 4 mm pockets for the magnet with its coil wound round it
  (glued as one: the model has only 0.5 mm between magnet and coil, too thin for a printed wall),
  a wire passage under each, four 8 mm gap-sensor holes beside the magnets towards the middle (no
  room outboard; an outboard first try left a 0.25 mm wall and a non-manifold mesh), and four
  3.03 mm skids inside the board's width. Probes: pocket floor, pocket wall and skid top at
  0.000 mm.
- Hull: an open tray (2.5 mm walls) for battery and electronics under the deck.
- Hangers: a constant cross-section printed standing on end, so every wall is vertical: a top
  plate over the 60 x 3 mm steel bar, tabs under the stator board's edges, side lips 2 mm off the
  cart's sides, two mounting holes; the station version adds ledges where the cart rests at the
  6.5 mm lift-off gap. The mounting holes are diamonds (round holes running sideways would have
  flat tops); the two remaining flags are the diamonds' apexes, which the mesher blunts into
  strips about 0.4 mm wide (located at the hole tops), not a real overhang.
- Mass: deck and hull 140 g solid, about 56 g at a typical 40% effective density; magnets
  4 x 3.4 g; coil copper 4 x 2.9 g. With the 200 g nominal payload that leaves about 220 g of the
  model's 500 g for battery, electronics, sensors and wiring.

### Parts list (demonstrator, one cart, 1.4 m of track)

| Item | Qty | Spec | From |
|---|---|---|---|
| NdFeB disc magnet | 4 | N42, 12 mm x 4 mm, axially magnetized (Br ~1.3 T) | `stray` |
| Control coil | 4 | 80 turns AWG 28 enamelled copper, 8 layers x 10, wound on a 13.0 mm mandrel and bonded; ~4 m and 0.86 ohm each | `corner` |
| Gap sensor | 4 | 3-8 mm range, noise about 1 um at 2 kHz (an inductive / eddy-current displacement sensor class) | `corner` noise sweep |
| Sideways sensor | 2 | the same class, reading the distance to the guide lips at x = +-60 mm | `cart` |
| Coil driver | 4 | H-bridge, +-12 V, 10 A peak, current sensing for a 1 kHz current loop | `corner` |
| Controller | 1 | 2 kHz gap loop (modal), current loops, zero-power loops | `cart` |
| Battery | 1 | 12 V class; about 0.2-0.4 W average (coils) plus electronics, 10 A peaks | `cart` |
| Steel flat bar | 1.4 m | low-carbon, 60 mm x 3 mm or thicker | `stray` (0.55 mm carries the flux; the image model wants it thick) |
| Stator board | 1.4 m | layer 1: three-phase motor coils, pitch 100/7 mm, 14.3 x 24 mm, 20 turns, two tracks 50 mm apart with phases reversed; layer 2: three-phase guidance coils 14.3 x 50 mm, 20 turns, spanning both rows; 1.5 mm thick | `motor`, `guide` |
| Motor driver | 1 | three-phase, block-switched sections under the cart, about 1 A per phase, position feedback for commutation | `motor`, `cart` |
| Guidance driver | 1 | three-phase, two independently driven blocks under the cart (front, rear), 8 A peak | `guide`, `cart` |
| Printed parts | | `cart_deck`, `cart_hull`, 14 x `rail_hanger`, 1-2 x `station_hanger` | `parts` |
| Skid pads | 4 | low-friction (PTFE) on the 3.03 mm skids | `parts` |

### Open items

- Sideways guidance is handled in simulation by the guidance coil layer (above); the side lips
  stay as a backstop for shoves well beyond 0.3 N at the cart's end.
- Both board layers were modelled with AWG 24-equivalent copper (0.13 and 0.22 ohm a coil);
  printed copper will be thinner, so recompute the losses for a real board.
- Power to the cart: a battery today; the wireless-charging project could feed it through the
  track.
- Hospital materials: PETG or polycarbonate and sealed seams for cleaning; not yet chosen.

## Plan

1. Forces, sizing and stray field (done for the demonstrator).
2. One corner in CortenForge's simulator (done): instability at the predicted rate, gap control
   through the coil, zero-power control, the gap-sensor requirement.
3. The cart (done): four corners, modal zero-power control, an off-centre payload, a stop at a
   station.
4. Linear-motor propulsion (done): a fine-pitch flat stator pushing on the lift magnets, in the
   cart trip.
5. Printable parts and a parts list (done). Open items listed above.
6. Sideways guidance (done): a guidance coil layer, front and rear blocks, sideways and yaw
   control in the six-motion cart.
