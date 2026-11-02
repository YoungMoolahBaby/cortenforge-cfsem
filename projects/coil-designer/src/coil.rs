//! Real coils: layer-wound round wire in a rectangular winding pack, one cfsem filament per turn.

use std::f64::consts::PI;

use cfsem::physics::circular_filament::{
    flux_circular_filament_scalar, flux_density_circular_filament_par,
    vector_potential_circular_filament_scalar,
};

/// (ohm m) annealed copper at 20 C, the IACS standard 1/58 uohm m.
pub const RHO_CU: f64 = 1.7241e-8;
/// (kg/m^3) copper.
pub const DENSITY_CU: f64 = 8960.0;
/// (m) assumed enamel build (both sides) on magnet wire; grade 2 is roughly this for AWG 18-26.
pub const ENAMEL_BUILD: f64 = 0.06e-3;

/// Enamelled round copper wire.
#[derive(Clone, Copy, Debug)]
pub struct Wire {
    pub awg: u32,
    /// (m) copper diameter.
    pub bare: f64,
    /// (m) diameter over enamel, which sets the winding pitch.
    pub outer: f64,
}

impl Wire {
    /// AWG by its definition, d = 0.127 mm * 92^((36 - n) / 39), plus [ENAMEL_BUILD].
    pub fn awg(n: u32) -> Self {
        let bare = 0.127e-3 * 92f64.powf((36.0 - n as f64) / 39.0);
        Wire {
            awg: n,
            bare,
            outer: bare + ENAMEL_BUILD,
        }
    }

    /// (ohm/m) at 20 C.
    pub fn ohm_per_m(&self) -> f64 {
        RHO_CU / (PI * self.bare * self.bare / 4.0)
    }
}

/// A rectangular winding pack, square-packed: `layers` radially, `per_layer` turns axially.
#[derive(Clone, Copy, Debug)]
pub struct Pack {
    /// (m) inner radius of the winding (the former's surface).
    pub r_in: f64,
    /// (m) axial centre.
    pub z_mid: f64,
    pub layers: usize,
    pub per_layer: usize,
    pub wire: Wire,
    /// +1 or -1: current direction relative to the series current.
    pub sign: f64,
}

impl Pack {
    pub fn turns(&self) -> usize {
        self.layers * self.per_layer
    }
    /// (m) radial depth.
    pub fn depth(&self) -> f64 {
        self.layers as f64 * self.wire.outer
    }
    /// (m) axial width.
    pub fn width(&self) -> f64 {
        self.per_layer as f64 * self.wire.outer
    }
    /// (m) mean winding radius.
    pub fn r_mean(&self) -> f64 {
        self.r_in + 0.5 * self.depth()
    }
    /// (r, z) of each turn's centre.
    pub fn centres(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        let p = self.wire.outer;
        let z0 = self.z_mid - 0.5 * self.width();
        (0..self.layers).flat_map(move |l| {
            (0..self.per_layer)
                .map(move |k| (self.r_in + (l as f64 + 0.5) * p, z0 + (k as f64 + 0.5) * p))
        })
    }
}

/// Series-connected packs.
#[derive(Clone, Debug)]
pub struct Coil {
    pub packs: Vec<Pack>,
}

impl Coil {
    /// One filament per turn: (r, z, sign, wire radius).
    pub fn filaments(&self) -> Vec<(f64, f64, f64, f64)> {
        self.packs
            .iter()
            .flat_map(|p| {
                p.centres()
                    .map(move |(r, z)| (r, z, p.sign, 0.5 * p.wire.bare))
            })
            .collect()
    }

    /// (Br, Bz) at the points for a series current (A), treating each turn as an ideal filament,
    /// which is exact outside the wire (a round wire's external field is its centreline's).
    pub fn field(&self, current: f64, r: &[f64], z: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let f = self.filaments();
        let rf: Vec<f64> = f.iter().map(|x| x.0).collect();
        let zf: Vec<f64> = f.iter().map(|x| x.1).collect();
        let nif: Vec<f64> = f.iter().map(|x| x.2 * current).collect();
        let wire = vec![0.0; f.len()];
        let (mut br, mut bz) = (vec![0.0; r.len()], vec![0.0; r.len()]);
        flux_density_circular_filament_par((&rf, &zf, &nif), &wire, (r, z), (&mut br, &mut bz))
            .expect("lengths match");
        (br, bz)
    }

    /// (m) total wire length (turns only, no leads).
    pub fn wire_length(&self) -> f64 {
        self.filaments().iter().map(|f| 2.0 * PI * f.0).sum()
    }

    /// (ohm) at 20 C.
    pub fn resistance(&self) -> f64 {
        self.packs
            .iter()
            .map(|p| p.centres().map(|(r, _)| 2.0 * PI * r).sum::<f64>() * p.wire.ohm_per_m())
            .sum()
    }

    /// (kg) copper.
    pub fn copper_mass(&self) -> f64 {
        self.packs
            .iter()
            .map(|p| {
                let len: f64 = p.centres().map(|(r, _)| 2.0 * PI * r).sum();
                len * PI * p.wire.bare * p.wire.bare / 4.0 * DENSITY_CU
            })
            .sum()
    }

    /// (H) series inductance, L = sum_ij s_i s_j M_ij. Mutual terms from cfsem's ideal-filament
    /// vector potential (exact between round wires, whose mutual inductance is that of their
    /// centrelines); self terms from cfsem's finite-wire model, see [ring_self_inductance].
    pub fn inductance(&self) -> f64 {
        let f = self.filaments();
        let mut cache: Vec<((f64, f64), f64)> = Vec::new();
        let mut l = 0.0;
        for (i, &(ri, zi, si, bi)) in f.iter().enumerate() {
            // Turns in a layer share a radius, so self terms repeat.
            let self_l = match cache.iter().find(|(k, _)| *k == (ri, bi)) {
                Some(&(_, v)) => v,
                None => {
                    let v = ring_self_inductance(ri, bi);
                    cache.push(((ri, bi), v));
                    v
                }
            };
            l += self_l;
            for &(rj, zj, sj, _) in &f[i + 1..] {
                let m = 2.0
                    * PI
                    * rj
                    * vector_potential_circular_filament_scalar((ri, zi, 1.0), (rj, zj));
                l += 2.0 * si * sj * m;
            }
        }
        l
    }
}

/// (H) self-inductance of one ring (radius a) of round wire (radius b), uniform current: the
/// mean over the wire's cross-section of cfsem's finite-wire flux per amp (energy method).
pub fn ring_self_inductance(a: f64, b: f64) -> f64 {
    ring_self_inductance_n(a, b, 16)
}

/// [ring_self_inductance] with `nt` angular points. Inside the wire cfsem's near-field potential
/// is a cubic polynomial in the local coordinates (physics/circular_filament.rs, the `t <= 1`
/// branch of `vector_potential_circular_filament_finite_thickness_scalar_near`), so the flux
/// 2 pi r A is a quartic: 4-point Gauss-Legendre in the radius (exact to degree 7, the area weight
/// included) and a uniform angular grid of more than 4 points integrate it exactly.
pub fn ring_self_inductance_n(a: f64, b: f64, nt: usize) -> f64 {
    // Gauss-Legendre on [0, 1]: nodes (1 +- x)/2, weights w/2.
    const GL4: [(f64, f64); 4] = [
        (0.069_431_844_202_973_71, 0.173_927_422_568_726_9),
        (0.330_009_478_207_571_9, 0.326_072_577_431_273_1),
        (0.669_990_521_792_428_1, 0.326_072_577_431_273_1),
        (0.930_568_155_797_026_3, 0.173_927_422_568_726_9),
    ];
    let mut sum = 0.0;
    for &(x, w) in &GL4 {
        let rho = x * b;
        for j in 0..nt {
            let th = (j as f64 + 0.5) / nt as f64 * 2.0 * PI;
            let p = (a + rho * th.cos(), rho * th.sin());
            sum += w * rho * flux_circular_filament_scalar((a, 0.0, 1.0), b, p);
        }
    }
    sum * b * (2.0 * PI / nt as f64) / (PI * b * b)
}

/// Maps (r, z) points to (Br, Bz).
pub type FieldFn<'a> = dyn Fn(&[f64], &[f64]) -> (Vec<f64>, Vec<f64>) + 'a;

/// Largest centred sphere radius (m) inside which |B - B0 z| / B0 <= eps, B0 the centre field.
/// The deviation is subharmonic in the bore, so its maximum over a ball sits on the surface,
/// and the surface max grows with the radius: bisection on the surface max is enough.
/// `field` maps (r, z) points to (Br, Bz).
pub fn uniform_radius(
    field: &FieldFn,
    eps: f64,
    r_max: f64,
) -> f64 {
    let b0 = field(&[0.0], &[0.0]).1[0];
    let n_theta = 91;
    let dev = |rho: f64| -> f64 {
        let (rs, zs): (Vec<f64>, Vec<f64>) = (0..n_theta)
            .map(|k| {
                let th = PI * k as f64 / (n_theta - 1) as f64;
                (rho * th.sin(), rho * th.cos())
            })
            .unzip();
        let (br, bz) = field(&rs, &zs);
        br.iter()
            .zip(&bz)
            .map(|(r, z)| r.hypot(z - b0) / b0.abs())
            .fold(0.0, f64::max)
    };
    let (mut lo, mut hi) = (0.0, r_max);
    if dev(hi) <= eps {
        return hi;
    }
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if dev(mid) <= eps { lo = mid } else { hi = mid }
    }
    lo
}
