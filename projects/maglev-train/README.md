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

## Plan

1. Forces, sizing and stray field (done for the demonstrator).
2. One corner in CortenForge's simulator (done): instability at the predicted rate, gap control
   through the coil, zero-power control, the gap-sensor requirement.
3. The cart: four corners, pitch and roll, payload changes, a stop at a station.
4. Linear-motor propulsion along the track.
5. Printable track and cart parts, a parts list.
