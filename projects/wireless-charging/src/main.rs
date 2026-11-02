//! Main 4: wireless power, first for the maglev cart (inductive power from the track).
//!
//! `check` tests cfsem's straight-filament inductance and field against textbook results; `link`
//! sizes the track loop and pickup; `effects` the rail, induced voltages and stray field.

use wireless_charging::{check, effects, former, layout, link, pad};

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "check".into());
    match mode.as_str() {
        "check" => {
            if !check::run() {
                std::process::exit(1);
            }
        }
        "link" => link::run(),
        "effects" => effects::run(),
        "pad" => pad::run(),
        "layout" => layout::run(),
        "former" => former::run(),
        m => {
            eprintln!("unknown mode {m}; modes: check, link, effects, pad, layout, former");
            std::process::exit(2);
        }
    }
}
