//! An active (magnetic) radiation shield for the mothership's torus habitat: a confined
//! toroidal field between two windings round the tube, cfsem's field, energy and forces, and
//! cosmic-ray transmission by backtracing.

mod check;
mod design;
mod dose;
mod scan;
mod shield;
mod size;
mod trace;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        "transmission" => scan::transmission(),
        "size" => size::run(),
        "failsafe" => size::failsafe(),
        other => {
            eprintln!("unknown mode {other:?}; modes: check, transmission, size, failsafe");
            std::process::exit(2);
        }
    }
}
