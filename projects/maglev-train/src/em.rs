//! Electromagnet model: air-core pancake coils under a steel rail, through cfsem.
//!
//! The rail's flat underside is the plane z = 0 and the coil sits below it. For an ideal soft-iron
//! rail (mu -> infinity, thick, wide, unsaturated) the field below the plane is that of the coil
//! plus its mirror image across the plane, carrying current in the same direction (so tangential
//! H vanishes on the surface). The rail's pull on the coil is then the force between the coil and
//! its image, and the coil's inductance gains their mutual inductance.

use std::f64::consts::PI;

use cfsem::physics::circular_filament::{
    flux_circular_filament, flux_circular_filament_scalar, flux_density_circular_filament,
    vector_potential_circular_filament_scalar,
};

/// (ohm m) annealed copper at 20 C.
pub const RHO_CU: f64 = 1.7241e-8;
/// (m) assumed enamel build on magnet wire.
pub const ENAMEL_BUILD: f64 = 0.06e-3;

#[derive(Clone, Copy, Debug)]
pub struct Wire {
    pub awg: u32,
    pub bare: f64,
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
}

/// A pancake coil, axis vertical, square-packed turns: `layers` radially from `r_in`,
/// `per_layer` axially. Its top face sits `gap` below the rail.
#[derive(Clone, Copy, Debug)]
pub struct Coil {
    pub r_in: f64,
    pub layers: usize,
    pub per_layer: usize,
    pub wire: Wire,
}

impl Coil {
    pub fn turns(&self) -> usize {
        self.layers * self.per_layer
    }
    pub fn r_out(&self) -> f64 {
        self.r_in + self.layers as f64 * self.wire.outer
    }
    pub fn height(&self) -> f64 {
        self.per_layer as f64 * self.wire.outer
    }
    /// (r, z) of each turn with the coil's top `gap` below the plane z = 0.
    pub fn turns_at(&self, gap: f64) -> Vec<(f64, f64)> {
        let p = self.wire.outer;
        let mut v = Vec::with_capacity(self.turns());
        for l in 0..self.layers {
            for k in 0..self.per_layer {
                v.push((
                    self.r_in + (l as f64 + 0.5) * p,
                    -gap - (k as f64 + 0.5) * p,
                ));
            }
        }
        v
    }
    /// (ohm) at 20 C.
    pub fn resistance(&self) -> f64 {
        let len: f64 = self.turns_at(0.0).iter().map(|t| 2.0 * PI * t.0).sum();
        RHO_CU * len / (PI * self.wire.bare * self.wire.bare / 4.0)
    }
}

/// An axially magnetized disc magnet as equivalent current loops: a uniform magnetization M
/// equals a surface current K = M = Br / mu0 around its side (NdFeB has mu_r ~ 1.05, so this is
/// close to exact). `n` loops spread over the height carry Br h / (mu0 n) each. Top face `gap`
/// below the plane. Returns (r, z) loops and the current per loop.
pub fn magnet_loops(
    radius: f64,
    height: f64,
    br: f64,
    n: usize,
    gap: f64,
) -> (Vec<(f64, f64)>, f64) {
    let loops = (0..n)
        .map(|k| (radius, -gap - (k as f64 + 0.5) * height / n as f64))
        .collect();
    (loops, br * height / (cfsem::MU_0 * n as f64))
}

/// Mirror images of turns across z = 0 (same current direction: an iron surface).
pub fn images(turns: &[(f64, f64)]) -> Vec<(f64, f64)> {
    turns.iter().map(|&(r, z)| (r, -z)).collect()
}

fn split(v: &[(f64, f64)]) -> (Vec<f64>, Vec<f64>) {
    v.iter().copied().unzip()
}

/// (N) upward force on coil `a` (current `ia` per turn) from loops `b` (current `ib` per turn),
/// by the Lorentz force on each turn: F_z = -2 pi r I B_r.
pub fn force_lorentz(a: &[(f64, f64)], ia: f64, b: &[(f64, f64)], ib: f64) -> f64 {
    let (rb, zb) = split(b);
    let cur = vec![ib; b.len()];
    let (ra, za) = split(a);
    let (mut br, mut bz) = (vec![0.0; a.len()], vec![0.0; a.len()]);
    flux_density_circular_filament(
        (&rb, &zb, &cur),
        &vec![0.0; b.len()],
        (&ra, &za),
        (&mut br, &mut bz),
    )
    .expect("lengths");
    ra.iter()
        .zip(&br)
        .map(|(r, b)| -2.0 * PI * r * ia * b)
        .sum()
}

/// (H) mutual inductance between two sets of series turns: flux of `b` (1 A per turn) linked by `a`.
pub fn mutual(a: &[(f64, f64)], b: &[(f64, f64)]) -> f64 {
    let (rb, zb) = split(b);
    let (ra, za) = split(a);
    let mut psi = vec![0.0; a.len()];
    flux_circular_filament(
        (&rb, &zb, &vec![1.0; b.len()]),
        &vec![0.0; b.len()],
        (&ra, &za),
        &mut psi,
    )
    .expect("lengths");
    psi.iter().sum()
}

/// (H) self-inductance of a coil in free space: turn self terms from cfsem's finite-wire flux
/// averaged over the wire section (exact 4-point Gauss quadrature for cfsem's interior model),
/// plus all mutual terms between turns from the ideal-filament potential.
pub fn self_inductance(c: &Coil) -> f64 {
    let t = c.turns_at(0.0);
    let b = 0.5 * c.wire.bare;
    let mut l = 0.0;
    let mut cache: Vec<(f64, f64)> = Vec::new();
    for (i, &(ri, zi)) in t.iter().enumerate() {
        let s = match cache.iter().find(|x| x.0 == ri) {
            Some(x) => x.1,
            None => {
                let v = ring_self(ri, b);
                cache.push((ri, v));
                v
            }
        };
        l += s;
        for &(rj, zj) in &t[i + 1..] {
            l += 2.0
                * 2.0
                * PI
                * rj
                * vector_potential_circular_filament_scalar((ri, zi, 1.0), (rj, zj));
        }
    }
    l
}

/// (H) one ring of round wire: the section mean of cfsem's finite-wire flux per amp.
fn ring_self(a: f64, b: f64) -> f64 {
    const GL4: [(f64, f64); 4] = [
        (0.069_431_844_202_973_71, 0.173_927_422_568_726_9),
        (0.330_009_478_207_571_9, 0.326_072_577_431_273_1),
        (0.669_990_521_792_428_1, 0.326_072_577_431_273_1),
        (0.930_568_155_797_026_3, 0.173_927_422_568_726_9),
    ];
    let nt = 16;
    let mut sum = 0.0;
    for &(x, w) in &GL4 {
        let rho = x * b;
        for j in 0..nt {
            let th = (j as f64 + 0.5) / nt as f64 * 2.0 * PI;
            sum += w
                * rho
                * flux_circular_filament_scalar(
                    (a, 0.0, 1.0),
                    b,
                    (a + rho * th.cos(), rho * th.sin()),
                );
        }
    }
    sum * b * (2.0 * PI / nt as f64) / (PI * b * b)
}

/// The rail's pull on a coil at `gap` with `current` (A), by the Lorentz force from the image.
pub fn lift(c: &Coil, gap: f64, current: f64) -> f64 {
    let t = c.turns_at(gap);
    force_lorentz(&t, current, &images(&t), current)
}

/// (H) the coil's mutual inductance with its image at `gap`; the coil's inductance near the rail
/// is L0 + this.
pub fn image_mutual(c: &Coil, gap: f64) -> f64 {
    let t = c.turns_at(gap);
    mutual(&t, &images(&t))
}
