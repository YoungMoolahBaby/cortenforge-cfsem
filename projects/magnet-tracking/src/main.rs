//! Main 1: track a permanent magnet's position and orientation from magnetometer readings.
//!
//! `check` tests the forward models (cfsem's dipole and current loops) against textbook results.

mod calib;
mod check;
mod layout;
mod magnet;
mod path;
mod solve;
mod tracker;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        "solve" => solve::run(),
        "layout" => layout::run(),
        "calib" => calib::run(),
        "path" => path::run(),
        m => {
            eprintln!("unknown mode {m}; modes: check, solve, layout, calib, path");
            std::process::exit(2);
        }
    }
}
