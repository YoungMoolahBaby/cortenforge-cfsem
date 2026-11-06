//! Free-space galactic cosmic-ray dose equivalent in tissue, and the share a shield lets
//! through. Spectra: the BPH00 local interstellar spectrum per nucleon (Usoskin et al. 2005,
//! eq. 2) modulated by the force-field approximation (eq. 1); element groups scaled by PDG's
//! abundances at 10.6 GeV/nucleon. Dose: Bethe stopping power in water (PDG eq. 34.5, no density
//! correction), ICRP 60 quality factor Q(L). Primaries only, in a thin target (skin-like): no
//! fragmentation or secondaries.

use std::f64::consts::PI;

/// Nucleon rest energy (GeV), electron rest energy (MeV), Bethe K (MeV mol^-1 cm^2), water's
/// mean excitation energy (MeV) and Z/A.
const TR: f64 = 0.938;
const ME: f64 = 0.510_998_95;
const K: f64 = 0.307_075;
const I_WATER: f64 = 79.7e-6;
const ZA_WATER: f64 = 0.555_09;
const AMU_MEV: f64 = 931.494;
const MEV_J: f64 = 1.602_176_634e-13;
pub const YEAR_S: f64 = 365.25 * 86400.0;

/// An element group: charge, mass number, abundance relative to oxygen at 10.6 GeV/nucleon
/// (PDG 2020 Table 30.1).
#[derive(Clone, Copy, Debug)]
pub struct Species {
    pub name: &'static str,
    pub z: f64,
    pub a: f64,
    pub abundance: f64,
}

pub const SPECIES: [Species; 12] = [
    Species { name: "H", z: 1.0, a: 1.0, abundance: 550.0 },
    Species { name: "He", z: 2.0, a: 4.0, abundance: 34.0 },
    Species { name: "Li-B", z: 4.0, a: 9.0, abundance: 0.40 },
    Species { name: "C-O", z: 7.0, a: 14.0, abundance: 2.20 },
    Species { name: "F-Ne", z: 10.0, a: 20.0, abundance: 0.30 },
    Species { name: "Na-Mg", z: 12.0, a: 24.0, abundance: 0.22 },
    Species { name: "Al-Si", z: 14.0, a: 28.0, abundance: 0.19 },
    Species { name: "P-S", z: 16.0, a: 32.0, abundance: 0.03 },
    Species { name: "Cl-Ar", z: 18.0, a: 37.0, abundance: 0.01 },
    Species { name: "K-Ca", z: 20.0, a: 40.0, abundance: 0.02 },
    Species { name: "Sc-Mn", z: 23.0, a: 50.0, abundance: 0.05 },
    Species { name: "Fe-Ni", z: 26.0, a: 56.0, abundance: 0.12 },
];

/// BPH00 local interstellar spectrum, particles / (m^2 sr s GeV/nucleon), T in GeV/nucleon.
pub fn lis(t: f64) -> f64 {
    let p = (t * (t + 2.0 * TR)).sqrt();
    1.9e4 * p.powf(-2.78) / (1.0 + 0.4866 * p.powf(-2.51))
}

/// Force-field modulated spectrum of a species at heliocentric modulation phi (GV).
pub fn spectrum(s: &Species, t: f64, phi_gv: f64) -> f64 {
    let big_phi = s.z / s.a * phi_gv;
    let scale = s.abundance / SPECIES[0].abundance;
    scale * lis(t + big_phi) * t * (t + 2.0 * TR) / ((t + big_phi) * (t + big_phi + 2.0 * TR))
}

/// Rigidity (GV) of a nucleus of the species at kinetic energy t (GeV/nucleon).
pub fn rigidity(s: &Species, t: f64) -> f64 {
    s.a / s.z * (t * (t + 2.0 * TR)).sqrt()
}

/// Mass stopping power in water (MeV cm^2 / g), Bethe without the density correction.
pub fn stopping(s: &Species, t: f64) -> f64 {
    let gamma = 1.0 + t / TR;
    let b2 = 1.0 - 1.0 / (gamma * gamma);
    let bg2 = b2 * gamma * gamma;
    let m = s.a * AMU_MEV;
    let wmax = 2.0 * ME * bg2 / (1.0 + 2.0 * gamma * ME / m + (ME / m).powi(2));
    K * s.z * s.z * ZA_WATER / b2 * (0.5 * (2.0 * ME * bg2 * wmax / (I_WATER * I_WATER)).ln() - b2)
}

/// ICRP 60 quality factor at unrestricted LET L (keV/um in water).
pub fn quality(l: f64) -> f64 {
    if l < 10.0 {
        1.0
    } else if l <= 100.0 {
        0.32 * l - 2.2
    } else {
        300.0 / l.sqrt()
    }
}

/// Free-space dose equivalent rate (Sv/yr) of a species, weighted by a transmission function of
/// rigidity (GV), and the absorbed dose (Gy/yr).
pub fn species_dose(s: &Species, phi_gv: f64, transmit: &dyn Fn(f64) -> f64) -> (f64, f64) {
    let (lo, hi, n) = (0.01f64.ln(), 1.0e4f64.ln(), 4000);
    let h = (hi - lo) / n as f64;
    let (mut dose, mut equiv) = (0.0, 0.0);
    for i in 0..=n {
        let t = (lo + h * i as f64).exp();
        let w = if i == 0 || i == n { 0.5 } else { 1.0 };
        let flux = 4.0 * PI * spectrum(s, t, phi_gv) * transmit(rigidity(s, t)); // m^-2 s^-1 (GeV/n)^-1
        let st = stopping(s, t); // MeV cm^2 / g
        let gy = flux * st * MEV_J * 1e3 * 1e-4 * t * h * w; // Gy/s over this d(ln T)
        dose += gy;
        equiv += gy * quality(0.1 * st);
    }
    (equiv * YEAR_S, dose * YEAR_S)
}

/// Total over the species: (Sv/yr, Gy/yr).
pub fn total(phi_gv: f64, transmit: &dyn Fn(f64) -> f64) -> (f64, f64) {
    SPECIES.iter().map(|s| species_dose(s, phi_gv, transmit)).fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1))
}

/// A transmission function from a table of (GV, fraction), linear in log rigidity, 0 below the
/// table and 1 above it.
pub fn from_table(table: Vec<(f64, f64)>) -> impl Fn(f64) -> f64 {
    move |r: f64| {
        if r <= table[0].0 {
            return table[0].1 * r / table[0].0;
        }
        for w in table.windows(2) {
            if r <= w[1].0 {
                let f = (r.ln() - w[0].0.ln()) / (w[1].0.ln() - w[0].0.ln());
                return w[0].1 + f * (w[1].1 - w[0].1);
            }
        }
        let last = table[table.len() - 1];
        last.1 + (1.0 - last.1) * (1.0 - last.0 / r)
    }
}
