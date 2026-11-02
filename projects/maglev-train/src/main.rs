//! Main 2: electromagnetic-suspension (EMS) maglev.
//!
//! `check` tests the force and inductance model (cfsem coils, iron rail by images).

mod check;
mod em;

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
