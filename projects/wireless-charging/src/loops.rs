//! Closed polyline coils and their inductances through cfsem's linear filaments.

use std::f64::consts::PI;

use cfsem::MU_0;
use cfsem::physics::linear_filament::{
    flux_density_linear_filament, inductance_piecewise_linear_filaments,
};

pub type P3 = (f64, f64, f64);

/// A closed polyline loop (last point joins the first), lumped as `turns` turns.
#[derive(Clone, Debug)]
pub struct Loop {
    pub pts: Vec<P3>,
    pub turns: f64,
}

impl Loop {
    /// A rectangle in the plane z, centred at (cx, cy), sides lx (along x) and ly, each side split
    /// into `n` segments, counter-clockwise seen from +z.
    pub fn rect(cx: f64, cy: f64, z: f64, lx: f64, ly: f64, n: usize, turns: f64) -> Self {
        let c = [
            (cx - lx / 2.0, cy - ly / 2.0),
            (cx + lx / 2.0, cy - ly / 2.0),
            (cx + lx / 2.0, cy + ly / 2.0),
            (cx - lx / 2.0, cy + ly / 2.0),
        ];
        let mut pts = Vec::new();
        for k in 0..4 {
            let (a, b) = (c[k], c[(k + 1) % 4]);
            for s in 0..n {
                let t = s as f64 / n as f64;
                pts.push((a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1), z));
            }
        }
        Loop { pts, turns }
    }

    /// A rectangle like `rect`, each side split into pieces no longer than `piece`.
    pub fn rect_by(cx: f64, cy: f64, z: f64, lx: f64, ly: f64, piece: f64, turns: f64) -> Self {
        let c = [
            (cx - lx / 2.0, cy - ly / 2.0),
            (cx + lx / 2.0, cy - ly / 2.0),
            (cx + lx / 2.0, cy + ly / 2.0),
            (cx - lx / 2.0, cy + ly / 2.0),
        ];
        let mut pts = Vec::new();
        for k in 0..4 {
            let (a, b) = (c[k], c[(k + 1) % 4]);
            let n = ((b.0 - a.0).hypot(b.1 - a.1) / piece).ceil().max(1.0) as usize;
            for s in 0..n {
                let t = s as f64 / n as f64;
                pts.push((a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1), z));
            }
        }
        Loop { pts, turns }
    }

    /// A circle of radius r centred at c, axis along +z, as an m-gon.
    pub fn circle(c: P3, r: f64, m: usize, turns: f64) -> Self {
        let pts = (0..m)
            .map(|k| {
                let t = 2.0 * PI * k as f64 / m as f64;
                (c.0 + r * t.cos(), c.1 + r * t.sin(), c.2)
            })
            .collect();
        Loop { pts, turns }
    }

    /// The mirror image across the plane z = 0: points reflected, order kept, so horizontal
    /// currents run the same way. Ideal soft iron: drive it with the same current; a perfect
    /// conductor: with the opposite current.
    pub fn image(&self) -> Self {
        Loop {
            pts: self.pts.iter().map(|p| (p.0, p.1, -p.2)).collect(),
            turns: self.turns,
        }
    }

    /// The mirror image across the plane z = p (order kept, as `image`).
    pub fn mirror(&self, p: f64) -> Self {
        Loop {
            pts: self.pts.iter().map(|q| (q.0, q.1, 2.0 * p - q.2)).collect(),
            turns: self.turns,
        }
    }

    /// Moved by (dx, dy, dz) after turning by `tilt` about the y axis through (0, 0, 0).
    pub fn placed(&self, tilt: f64, d: P3) -> Self {
        let (s, c) = tilt.sin_cos();
        Loop {
            pts: self
                .pts
                .iter()
                .map(|q| (c * q.0 + s * q.2 + d.0, q.1 + d.1, -s * q.0 + c * q.2 + d.2))
                .collect(),
            turns: self.turns,
        }
    }

    /// (start, delta) segments.
    pub fn segments(&self) -> (Vec<P3>, Vec<P3>) {
        let n = self.pts.len();
        let deltas = (0..n)
            .map(|k| {
                let (a, b) = (self.pts[k], self.pts[(k + 1) % n]);
                (b.0 - a.0, b.1 - a.1, b.2 - a.2)
            })
            .collect();
        (self.pts.clone(), deltas)
    }

    /// (m) length of one turn.
    pub fn length(&self) -> f64 {
        let (_, d) = self.segments();
        d.iter()
            .map(|v| (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt())
            .sum()
    }
}

fn split3(v: &[P3]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    (
        v.iter().map(|p| p.0).collect(),
        v.iter().map(|p| p.1).collect(),
        v.iter().map(|p| p.2).collect(),
    )
}

/// (H) mutual inductance between two loops (turns included), cfsem's A.dl kernel.
pub fn mutual(a: &Loop, b: &Loop) -> f64 {
    let (sa, da) = a.segments();
    let (sb, db) = b.segments();
    let (sa, da, sb, db) = (split3(&sa), split3(&da), split3(&sb), split3(&db));
    let r = vec![1e-6; a.pts.len()];
    a.turns
        * b.turns
        * inductance_piecewise_linear_filaments(
            (&sa.0, &sa.1, &sa.2),
            (&da.0, &da.1, &da.2),
            (&sb.0, &sb.1, &sb.2),
            (&db.0, &db.1, &db.2),
            &r,
        )
        .expect("lengths")
}

/// (H) cfsem's self-inductance of one turn of round wire (radius `wire`), as cfsem returns it: its
/// A.dl runs on the wire axis, which counts the internal inductance as mu0/4pi per metre.
pub fn self_inductance_raw(a: &Loop, wire: f64) -> f64 {
    let (s, d) = a.segments();
    let (s, d) = (split3(&s), split3(&d));
    let r = vec![wire; a.pts.len()];
    inductance_piecewise_linear_filaments(
        (&s.0, &s.1, &s.2),
        (&d.0, &d.1, &d.2),
        (&s.0, &s.1, &s.2),
        (&d.0, &d.1, &d.2),
        &r,
    )
    .expect("lengths")
}

/// (H) self-inductance of the loop, turns squared, uniform (DC or Litz) current: cfsem's value less
/// mu0/8pi per metre, so the internal term is the cross-section average, mu0/8pi per metre.
pub fn self_inductance(a: &Loop, wire: f64) -> f64 {
    a.turns * a.turns * (self_inductance_raw(a, wire) - MU_0 / (8.0 * PI) * a.length())
}

/// (T) field of loops carrying currents (amps per turn) at points.
pub fn field(loops: &[(&Loop, f64)], pts: &[P3]) -> Vec<P3> {
    let (mut s, mut d, mut i) = (Vec::new(), Vec::new(), Vec::new());
    for (l, cur) in loops {
        let (ls, ld) = l.segments();
        s.extend(ls);
        d.extend(ld);
        i.extend(std::iter::repeat_n(cur * l.turns, l.pts.len()));
    }
    let (s, d, p) = (split3(&s), split3(&d), split3(pts));
    let zero = vec![0.0; i.len()];
    let (mut bx, mut by, mut bz) = (
        vec![0.0; pts.len()],
        vec![0.0; pts.len()],
        vec![0.0; pts.len()],
    );
    flux_density_linear_filament(
        (&p.0, &p.1, &p.2),
        (&s.0, &s.1, &s.2),
        (&d.0, &d.1, &d.2),
        &i,
        &zero,
        (&mut bx, &mut by, &mut bz),
    )
    .expect("lengths");
    (0..pts.len()).map(|k| (bx[k], by[k], bz[k])).collect()
}
