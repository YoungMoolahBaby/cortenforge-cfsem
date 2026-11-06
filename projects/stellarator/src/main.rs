//! A toy modular stellarator: Fourier-series coils, cfsem's Biot-Savart field, field lines and
//! Poincare sections, and the coils swept into printable solids with cf-design.

mod check;
mod coils;
mod design;
mod meshes;
mod poincare;
mod render;
mod scan;
mod trace;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        "scan" => scan::run(),
        "trace" => poincare::run(),
        "mesh" => meshes::run(),
        "mesh-probe" => meshes::probe(),
        "sdf-probe" => meshes::sdf_probe(),
        other => {
            eprintln!("unknown mode {other:?}; modes: check, scan, trace, mesh");
            std::process::exit(2);
        }
    }
}
