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

## Plan

1. Checks (done).
2. Track loop and pickup for the cart: coupling against sideways and vertical motion, the 85 kHz
   series-tuned link (power, efficiency, track current), the steel rail bounded by images.
3. Side effects: 85 kHz voltage induced in the lift coils (gap-control noise), stray field.
4. Printable pickup former and track tray; the pickup on the cart in the maglev simulation.
