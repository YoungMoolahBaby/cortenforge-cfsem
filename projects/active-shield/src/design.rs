//! The baseline shield for the mothership (C:\dev\mothership): a 70 m habitat tube whose centre
//! circle is 824 m from the spin axis, windings at 75 and 95 m, 5 T at the tube centre's radius
//! (100 T m of bending power).

use crate::shield::Shield;

/// Tube centre's distance from the spin axis (m): the mothership's 894 m floor less its 70 m tube.
pub const R_MAJOR: f64 = 824.0;
/// Habitat outer radius (m), the mothership's tube.
pub const R_HAB: f64 = 70.0;
/// Inner winding (m): clear of the habitat, its 1.5 m gap and a thin residual passive shield.
pub const R1: f64 = 75.0;
/// Shell thickness (m) and field (T) for the baseline.
pub const THICKNESS: f64 = 20.0;
pub const B_CENTRE: f64 = 5.0;
/// Loops a winding: ~4.7 m apart on the inner equator.
pub const LOOPS: usize = 1000;

pub fn baseline() -> Shield {
    Shield::with_field(R_MAJOR, R1, R1 + THICKNESS, LOOPS, B_CENTRE)
}
