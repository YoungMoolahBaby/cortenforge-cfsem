//! Warm-up 1: coil designer.
//!
//! `check` tests cfsem's circular-filament API against known results.

mod check;
use coil_designer::coil;
mod design;
mod optimize;
mod pack;
mod refs;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        "pack" => pack::run(),
        "optimize" => optimize::run(),
        "design" => design::run(),
        "former" => design::run_former(&std::env::args().skip(2).collect::<Vec<_>>()),
        other => {
            eprintln!("unknown mode {other:?}; modes: check, pack, optimize, design, former");
            std::process::exit(2);
        }
    }
}
