# cortenforge-cfsem

Example projects that join [CortenForge](https://crates.io/crates/cortenforge) (simulation, geometry,
optimization, meshes, fabrication) with [cfsem](https://crates.io/crates/cfsem) (quasi-steady
electromagnetics: filaments, Biot-Savart, inductance, Grad-Shafranov).

One Cargo workspace; each project lives in `projects/<name>/`.

## Projects

**Warm-ups**

| Project | What |
|---|---|
| [coil-designer](projects/coil-designer) | Optimize Helmholtz or Maxwell coils for a uniform field; the gentlest way to learn cfsem |
| [fusion-visualizer](projects/fusion-visualizer) | Field lines around a small tokamak or stellarator, and the plasma equilibrium shape |

**Main projects**

| Project | What |
|---|---|
| [magnet-tracking](projects/magnet-tracking) | Track a magnet from magnetometer boards; becomes the localization module for colon-probe |
| [maglev-train](projects/maglev-train) | Electromagnetic-suspension train on a printable tabletop track |
| [colon-probe](projects/colon-probe) | Magnetically steered colonoscopy probe: tethered tip, then an untethered capsule; Endotics-style inchworm benchmark |
| [wireless-charging](projects/wireless-charging) | Map the pad-to-receiver sweet spot and optimize the pad layout; could power the capsule |
| [soft-robots](projects/soft-robots) | Silicone robots with embedded magnets that bend, crawl or grip; cast with CortenForge mold tooling |
| [stellarator](projects/stellarator) | A tabletop vacuum-field stellarator: Fourier coils, nested surfaces from field lines, printable coil sweeps; coil search next |

**Extras**

| Project | What |
|---|---|
| [planar-conveyor](projects/planar-conveyor) | A coil grid that levitates and routes several movers |
| [microswimmers](projects/microswimmers) | Helical swimmers driven by rotating fields (depends on CortenForge fluid support) |
| [magnet-arrays](projects/magnet-arrays) | Halbach arrays and contactless magnetic gears, magnets as equivalent current loops |

## Build

```sh
cargo run --release -p coil-designer
```

Rust 2024 edition; cortenforge 0.9, cfsem 14.
