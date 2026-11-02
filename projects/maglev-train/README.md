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

- `check` (default): force and inductance model against known results, 16 checks, plus design
  numbers for air-core vs hybrid. Exits 1 on any failure.

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

## Plan

1. Forces (in progress): hybrid force vs gap and current, linearized stiffness, inductance, the
   stray-field map, sizing for the demonstrator.
2. One corner in CortenForge's simulator: the open-loop instability at its predicted rate, then
   gap control through the coil (voltage drive, inductance, sensor noise), then zero-power
   control.
3. The cart: four corners, pitch and roll, payload changes, a stop at a station.
4. Linear-motor propulsion along the track.
5. Printable track and cart parts, a parts list.
