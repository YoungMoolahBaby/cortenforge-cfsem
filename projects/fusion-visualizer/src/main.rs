//! Warm-up 2: fusion magnet visualizer.
//!
//! `check` tests cfsem's Grad-Shafranov operator and its flux and field conventions.

mod check;
mod gs;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        other => {
            eprintln!("unknown mode {other:?}; modes: check");
            std::process::exit(2);
        }
    }
}
