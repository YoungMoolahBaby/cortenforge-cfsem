//! Main 4: wireless power, first for the maglev cart (inductive power from the track).
//!
//! `check` tests cfsem's straight-filament inductance and field against textbook results.

mod check;
mod link;
mod loops;

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        "link" => link::run(),
        m => {
            eprintln!("unknown mode {m}; modes: check, link");
            std::process::exit(2);
        }
    }
}
