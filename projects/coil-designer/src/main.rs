//! Warm-up 1: coil designer.
//!
//! `check` tests cfsem's circular-filament API against known results.

mod check;

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
