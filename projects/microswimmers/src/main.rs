//! Extra: helical microswimmers turned by a rotating magnetic field.
//!
//! `check` tests the physics (regularized Stokeslets, resistive force theory, CortenForge's fluid
//! model, cfsem's coil fields) before anything is built on it.

mod check;
mod coils;
mod helix;
mod rft;
mod simcore;
mod stokeslet;
mod swim;
mod swimmer;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "swim" => swim::run(),
        "coils" => coils::run(),
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        m => {
            eprintln!("unknown mode {m}; modes: check, swim, coils");
            std::process::exit(2);
        }
    }
}
