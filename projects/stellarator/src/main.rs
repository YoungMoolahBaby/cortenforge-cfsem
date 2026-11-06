//! A toy modular stellarator: Fourier-series coils, cfsem's Biot-Savart field, field lines and
//! Poincare sections, and the coils swept into printable solids with cf-design.

mod check;
mod coils;
mod design;
mod meshes;
mod optimize;
mod poincare;
mod render;
mod scan;
mod search;
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
        "opt-probe" => search::probe(),
        "best" => {
            let rest: Vec<String> = std::env::args().skip(2).collect();
            optimize::best(&rest)
        }
        "optimize" => {
            let rest: Vec<String> = std::env::args().skip(2).collect();
            optimize::run(&rest)
        }
        "trace" => poincare::run(),
        "mesh" => meshes::run(),
        "mesh-probe" => meshes::probe(),
        "sdf-probe" => meshes::sdf_probe(),
        other => {
            eprintln!("unknown mode {other:?}; modes: check, scan, trace, mesh, opt-probe, optimize, best");
            std::process::exit(2);
        }
    }
}
