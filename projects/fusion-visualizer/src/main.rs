//! Warm-up 2: fusion magnet visualizer.
//!
//! `check` tests cfsem's Grad-Shafranov operator and its flux and field conventions.

mod check;
mod contour;
mod equilibrium;
mod gs;
mod plot;
mod shafranov;
mod shape;
mod tf;
mod trace;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        "shafranov" => shafranov::run(),
        "shape" => shape::run(),
        "trace" => trace::run(),
        other => {
            eprintln!("unknown mode {other:?}; modes: check, shafranov, shape, trace");
            std::process::exit(2);
        }
    }
}
