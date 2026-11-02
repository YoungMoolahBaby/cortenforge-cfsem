//! Contours of a grid field (marching squares) and helpers for closed curves.

use crate::gs::Grid;

/// Polylines where `f` crosses `level`. Segments are linked into polylines; closed ones repeat
/// their first point at the end.
pub fn contours(grid: &Grid, f: &[f64], level: f64) -> Vec<Vec<(f64, f64)>> {
    let (nr, nz) = (grid.nr(), grid.nz());
    let at = |i: usize, j: usize| f[grid.idx(i, j)] - level;
    let lerp = |p: (f64, f64), q: (f64, f64), a: f64, b: f64| {
        let t = a / (a - b);
        (p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1))
    };
    let mut segs: Vec<((f64, f64), (f64, f64))> = Vec::new();
    for i in 0..nr - 1 {
        for j in 0..nz - 1 {
            let c = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)];
            let v: Vec<f64> = c.iter().map(|&(a, b)| at(a, b)).collect();
            if v.iter().any(|x| !x.is_finite()) {
                continue;
            }
            let p: Vec<(f64, f64)> = c.iter().map(|&(a, b)| (grid.rs[a], grid.zs[b])).collect();
            // Crossing points on the four edges.
            let mut pts = Vec::new();
            for e in 0..4 {
                let (a, b) = (e, (e + 1) % 4);
                if (v[a] < 0.0) != (v[b] < 0.0) {
                    pts.push(lerp(p[a], p[b], v[a], v[b]));
                }
            }
            match pts.len() {
                2 => segs.push((pts[0], pts[1])),
                // Saddle cell: pair by the cell-centre value.
                4 => {
                    let centre = v.iter().sum::<f64>() / 4.0;
                    if (centre < 0.0) == (v[0] < 0.0) {
                        segs.push((pts[0], pts[3]));
                        segs.push((pts[1], pts[2]));
                    } else {
                        segs.push((pts[0], pts[1]));
                        segs.push((pts[2], pts[3]));
                    }
                }
                _ => {}
            }
        }
    }
    link(segs)
}

/// Join segments that share end points into polylines.
fn link(segs: Vec<((f64, f64), (f64, f64))>) -> Vec<Vec<(f64, f64)>> {
    use std::collections::HashMap;
    let key = |p: (f64, f64)| ((p.0 * 1e9).round() as i64, (p.1 * 1e9).round() as i64);
    let mut ends: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (k, s) in segs.iter().enumerate() {
        ends.entry(key(s.0)).or_default().push(k);
        ends.entry(key(s.1)).or_default().push(k);
    }
    let mut used = vec![false; segs.len()];
    let mut lines = Vec::new();
    for start in 0..segs.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut line = vec![segs[start].0, segs[start].1];
        // Grow forward from the tail, then backward from the head.
        for forward in [true, false] {
            loop {
                let tip = if forward {
                    *line.last().expect("non-empty")
                } else {
                    line[0]
                };
                let next = ends
                    .get(&key(tip))
                    .and_then(|ks| ks.iter().copied().find(|&k| !used[k]));
                let Some(k) = next else { break };
                used[k] = true;
                let s = segs[k];
                let other = if key(s.0) == key(tip) { s.1 } else { s.0 };
                if forward {
                    line.push(other);
                } else {
                    line.insert(0, other);
                }
            }
        }
        lines.push(line);
    }
    lines
}

/// Even-odd point-in-polygon test.
pub fn inside(poly: &[(f64, f64)], p: (f64, f64)) -> bool {
    let mut c = false;
    let n = poly.len();
    for k in 0..n {
        let (a, b) = (poly[k], poly[(k + n - 1) % n]);
        if (a.1 > p.1) != (b.1 > p.1) && p.0 < (b.0 - a.0) * (p.1 - a.1) / (b.1 - a.1) + a.0 {
            c = !c;
        }
    }
    c
}

/// Perimeter of a closed polyline.
pub fn perimeter(poly: &[(f64, f64)]) -> f64 {
    poly.windows(2)
        .map(|w| (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1))
        .sum()
}
