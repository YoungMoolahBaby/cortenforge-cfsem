//! Extra: helical microswimmers turned by a rotating magnetic field.
//!
//! `check` tests the physics (regularized Stokeslets, resistive force theory, CortenForge's fluid
//! model, cfsem's coil fields) before anything is built on it.

mod check;
mod coils;
mod helix;
mod rft;
mod scale;
mod shape;
mod simcore;
mod steer;
mod stokeslet;
mod swim;
mod swimmer;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "swim" => swim::run(),
        "coils" => coils::run(),
        "scale" => scale::run(),
        "steer" => match std::env::args().nth(2).as_deref() {
            Some("wall") => steer::run_wall(),
            Some("contact") => steer::run_contact(),
            _ => steer::run(),
        },
        "shape" => match std::env::args().nth(2).as_deref() {
            Some("part") => shape::run_part(),
            Some("probe") => shape::run_probe(),
            _ => shape::run(),
        },
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        m => {
            eprintln!("unknown mode {m}; modes: check, swim, coils, shape, steer, scale");
            std::process::exit(2);
        }
    }
}
