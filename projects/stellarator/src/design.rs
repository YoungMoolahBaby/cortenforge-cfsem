//! The toy stellarator's numbers, shared by every mode.

use crate::coils::{Filaments, RotatingEllipse};
use crate::trace::Tracer;

/// (T) toroidally averaged field at R0.
pub const B0: f64 = 0.1;
/// (m) major radius.
pub const R0: f64 = 0.25;
/// Straight filaments per coil.
pub const SEGS: usize = 128;
/// (m) conductor radius (cfsem's near-wire regularization, and the swept pipe's radius).
pub const WIRE_RADIUS: f64 = 0.004;

pub fn baseline() -> RotatingEllipse {
    RotatingEllipse {
        n_fp: 3,
        coils_per_period: 6,
        rc: R0,
        a_long: 0.10,
        a_short: 0.04,
        shift: 0.02,
        shift_z: -0.02,
        current: 0.0,
    }
    .with_field(B0, R0)
}

pub fn filaments(c: &RotatingEllipse) -> Filaments {
    Filaments::new(&c.coils(), c.current, SEGS, WIRE_RADIUS)
}

pub fn tracer(c: &RotatingEllipse) -> Tracer {
    Tracer {
        n_fp: c.n_fp,
        steps: 120,
        centre: (c.rc, 0.0),
        lost_radius: c.a_long + c.shift,
    }
}
