//! Main 2: electromagnetic-suspension (EMS) maglev.
//!
//! `check` tests the force and inductance model (cfsem coils, iron rail by images).

mod cart;
mod check;
mod corner;
mod em;
mod guide;
mod motor;
mod parts;
mod stray;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        "stray" => stray::run(),
        "corner" => corner::run(),
        "cart" => cart::run(),
        "parts" => parts::run(),
        "guide" => {
            guide::run();
        }
        "motor" => {
            motor::run();
        }
        other => {
            eprintln!(
                "unknown mode {other:?}; modes: check, stray, corner, cart, motor, parts, guide"
            );
            std::process::exit(2);
        }
    }
}
