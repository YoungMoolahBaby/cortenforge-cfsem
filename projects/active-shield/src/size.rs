//! `size`: what each shield weighs and lets through. The conductor carries N I round both
//! windings at an engineering current density; the structure that holds the magnetic pressure
//! weighs Q rho E / sigma (the virial theorem gives Q >= 1); the inner winding, facing the warm
//! habitat, leaks heat that cryocoolers must lift. The dose through the magnet is the
//! free-space dose of what it transmits (those particles are all above ~10 GV and cross the body
//! whole, so skin and organ doses match). A hybrid adds water behind the magnet until the dose
//! meets the target, at the e-fold measured deep in water where the spectrum is similarly hard.

use std::f64::consts::PI;

use crate::design;
use crate::dose;
use crate::scan::{curve, rigidity_at, RIGIDITIES_GV};
use crate::shield::Shield;
use crate::trace::Tracer;

/// Inputs (README "Inputs" gives their sources).
pub struct Inputs {
    /// Engineering current density of the conductor (A/m^2) and its density (kg/m^3).
    pub j_e: f64,
    pub rho_conductor: f64,
    /// Structure: allowable stress (Pa), density (kg/m^3), virial factor Q.
    pub sigma: f64,
    pub rho_structure: f64,
    pub q: f64,
    /// Thin passive shield always kept inside the magnet (t/m^2), and the area it covers (m^2).
    pub residual: f64,
    pub shield_area: f64,
    /// Heat leak onto the inner winding (W/m^2) and cryocooler input per watt lifted at ~40 K.
    pub heat_leak: f64,
    pub cooler_w_per_w: f64,
    /// Secondaries the transmitted (GeV) primaries make in the hull, water and body, as a factor
    /// on the primaries-only dose (judgement).
    pub buildup: f64,
    /// Solar modulation (GV), the dose target (Sv/yr), and water's e-fold of dose (t/m^2).
    pub phi_gv: f64,
    pub target: f64,
    pub water_efold: f64,
}

pub fn inputs() -> Inputs {
    Inputs {
        j_e: 1.0e9,
        rho_conductor: 8870.0,
        sigma: 1.2e9,
        rho_structure: 1550.0,
        q: 3.0,
        residual: 0.5,
        // The mothership's still shield round its 70 m tube with a 1.5 m gap: 4 pi^2 R (r + gap).
        shield_area: 4.0 * PI * PI * design::R_MAJOR * (design::R_HAB + 1.5),
        heat_leak: 1.25,
        cooler_w_per_w: 43.0,
        buildup: 2.0,
        phi_gv: 0.42,
        target: 0.010,
        water_efold: 2.01,
    }
}

pub struct Mass {
    pub energy: f64,
    pub conductor: f64,
    pub structure: f64,
    pub cryo_w: f64,
}

impl Mass {
    pub fn magnet(&self) -> f64 {
        self.conductor + self.structure
    }
}

pub fn mass(sh: &Shield, inp: &Inputs) -> Mass {
    let amp_metres = sh.amp_turns * 2.0 * PI * (sh.r1 + sh.r2);
    let energy = sh.energy_ideal();
    let inner_area = 4.0 * PI * PI * sh.r_major * sh.r1;
    Mass {
        energy,
        conductor: amp_metres / inp.j_e * inp.rho_conductor,
        structure: inp.q * inp.rho_structure * energy / inp.sigma,
        cryo_w: inner_area * inp.heat_leak * inp.cooler_w_per_w,
    }
}

pub fn run() {
    let inp = inputs();
    let tr = Tracer { r_hab: design::R_HAB, ..Tracer::default() };
    let one = |_: f64| 1.0;
    let (free, _) = dose::total(inp.phi_gv, &one);
    println!(
        "Inputs: J_e {:.0} A/mm^2 at {:.0} kg/m^3; structure {:.0} MPa at {:.0} kg/m^3, Q = {:.1}; residual water {:.1} t/m^2 over {:.3e} m^2;",
        inp.j_e / 1e6, inp.rho_conductor, inp.sigma / 1e6, inp.rho_structure, inp.q, inp.residual, inp.shield_area
    );
    println!(
        "        heat leak {:.2} W/m^2 at {:.0} W/W; secondaries x{:.1}; solar minimum (phi {:.2} GV): free-space dose equivalent {:.0} mSv/yr (skin); target {:.0} mSv/yr\n",
        inp.heat_leak, inp.cooler_w_per_w, inp.buildup, inp.phi_gv, free * 1e3, inp.target * 1e3
    );
    println!(
        "  {:>4} {:>3} {:>5} | {:>5} | {:>8} {:>6} {:>6} {:>6} | {:>6} {:>6} | {:>6} {:>7} {:>7} | {:>5}",
        "T m", "d", "B T", "50%GV", "energy J", "cond", "struct", "magnet", "prim", "x2", "+water", "hybrid", "t/m^2", "cryoMW"
    );
    for (bl, d) in [(50.0, 20.0), (100.0, 20.0), (100.0, 40.0), (150.0, 40.0), (200.0, 40.0), (200.0, 60.0), (300.0, 60.0), (400.0, 80.0)] {
        let s = Shield::with_field(design::R_MAJOR, design::R1, design::R1 + d, design::LOOPS, bl / d);
        let c = curve(&s, &tr, &RIGIDITIES_GV);
        let table: Vec<(f64, f64)> = c.iter().map(|r| (r.0, r.1)).collect();
        let transmit = dose::from_table(table);
        let (through, _) = dose::total(inp.phi_gv, &transmit);
        let m = mass(&s, &inp);
        // Water behind the magnet to reach the target (at least the residual layer).
        let water = (inp.water_efold * (inp.buildup * through / inp.target).ln()).max(inp.residual);
        let water_t = water * 1e3 * inp.shield_area;
        let hybrid = m.magnet() + water_t;
        println!(
            "  {bl:>4.0} {d:>3.0} {:>5.2} | {:>5.1} | {:>8.2e} {:>6.2} {:>6.2} {:>6.2} | {:>6.0} {:>6.0} | {:>6.2} {:>7.2} {:>7.2} | {:>5.0}",
            bl / d,
            rigidity_at(&c, 0.5),
            m.energy,
            m.conductor / 1e9,
            m.structure / 1e9,
            m.magnet() / 1e9,
            through * 1e3,
            inp.buildup * through * 1e3,
            water,
            hybrid / 1e9,
            hybrid / 1e3 / inp.shield_area,
            m.cryo_w / 1e6
        );
    }
    println!(
        "\n  masses in Mt; doses mSv/yr through the magnet, primaries only and x{:.1} for secondaries; +water = t/m^2 behind it for {:.0} mSv/yr",
        inp.buildup,
        inp.target * 1e3
    );
    println!("  water alone (the mothership's design): 8.0 t/m^2 = {:.1} Mt", 8.0 * inp.shield_area / 1e6);
}


/// Water (t/m^2) for a yearly dose `limit` (Sv) with no magnet: the OLTARIS water curve of Globus
/// & Strout 2017 from its deepest row (7 t/m^2 gives 16.1 mSv/yr) at 2.01 t/m^2 an e-fold, the
/// mothership's own water model.
pub fn water_for(limit: f64) -> f64 {
    7.0 + 2.01 * (0.0161 / limit).ln()
}

/// `failsafe`: water thick enough that the dose stays under a limit with the magnet off, and the
/// smallest magnet that brings it to the target with the magnet on. Shells 40 m thick (lighter
/// and sharper for the same bending power). With the magnet on, the dose behind the water is the
/// magnet's leftovers (x the secondaries factor) attenuated at water's e-fold, never more than
/// the magnet-off dose. Cryocoolers at ~40 K (J_e 1,000 A/mm^2, 43 W/W) or, for these low-field
/// magnets, warmer at ~65 K (J_e ~300 A/mm^2, 18 W/W).
pub fn failsafe() {
    let inp = inputs();
    let tr = Tracer { r_hab: design::R_HAB, ..Tracer::default() };
    let d = 40.0;
    let candidates = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 75.0, 100.0, 125.0, 150.0];
    let warm_je = 3.0e8;
    let warm_w_per_w = 18.0;
    println!(
        "Fail-safe hybrids: water for a magnet-off limit, the smallest magnet (d = {d:.0} m) for {:.0} mSv/yr with it on; secondaries x{:.1}\n",
        inp.target * 1e3,
        inp.buildup
    );
    // Each candidate magnet's leftovers (primaries, Sv/yr), mass and energy.
    let mags: Vec<(f64, f64, Mass, Mass)> = candidates
        .iter()
        .map(|&bl| {
            let s = Shield::with_field(design::R_MAJOR, design::R1, design::R1 + d, design::LOOPS, bl / d);
            let c = curve(&s, &tr, &RIGIDITIES_GV);
            let transmit = dose::from_table(c.iter().map(|r| (r.0, r.1)).collect());
            let (through, _) = dose::total(inp.phi_gv, &transmit);
            let cold = mass(&s, &inp);
            let warm = mass(&s, &Inputs { j_e: warm_je, cooler_w_per_w: warm_w_per_w, ..inputs() });
            println!("  magnet {bl:>4.0} T m ({:.2} T): leftovers {:>6.1} mSv/yr (primaries), magnet {:.3} Mt, {:.2e} J", bl / d, through * 1e3, cold.magnet() / 1e9, cold.energy);
            (bl, through, cold, warm)
        })
        .collect();
    println!(
        "\n  {:>7} | {:>6} {:>6} | {:>5} {:>7} | {:>7} {:>7} | {:>8} {:>6} | {:>6} {:>6}",
        "off mSv", "water", "Mt", "T m", "magnet", "total", "on mSv", "energy J", "kt TNT", "40K MW", "65K MW"
    );
    for limit in [0.020, 0.050, 0.100] {
        let water = water_for(limit).max(inp.residual);
        let water_t = water * 1e3 * inp.shield_area;
        let on = |through: f64| (inp.buildup * through * (-water / inp.water_efold).exp()).min(limit);
        match mags.iter().find(|m| on(m.1) <= inp.target) {
            Some((bl, through, cold, warm)) => println!(
                "  {:>7.0} | {:>6.2} {:>6.2} | {:>5.0} {:>7.3} | {:>7.2} {:>7.1} | {:>8.2e} {:>6.0} | {:>6.0} {:>6.0}",
                limit * 1e3,
                water,
                water_t / 1e9,
                bl,
                cold.magnet() / 1e9,
                (water_t + cold.magnet()) / 1e9,
                on(*through) * 1e3,
                cold.energy,
                cold.energy / 4.184e12,
                cold.cryo_w / 1e6,
                warm.cryo_w / 1e6
            ),
            None => println!("  {:>7.0} | {:>6.2} {:>6.2} | no candidate magnet reaches the target", limit * 1e3, water, water_t / 1e9),
        }
    }
    println!(
        "\n  warm (65 K) magnets carry ~3x the conductor: e.g. {:.3} Mt instead of {:.3} Mt at 50 T m; water alone for {:.0} mSv/yr: {:.2} t/m^2 = {:.1} Mt",
        mags.iter().find(|m| m.0 == 50.0).map(|m| m.3.magnet() / 1e9).unwrap_or(f64::NAN),
        mags.iter().find(|m| m.0 == 50.0).map(|m| m.2.magnet() / 1e9).unwrap_or(f64::NAN),
        inp.target * 1e3,
        water_for(inp.target),
        water_for(inp.target) * inp.shield_area / 1e6
    );
}
