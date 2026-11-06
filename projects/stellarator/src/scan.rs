//! `scan`: axis and rotational transform across coil shapes, to pick the baseline.

use std::time::Instant;

use rayon::prelude::*;

use crate::coils::RotatingEllipse;
use crate::design;

pub fn run() {
    let base = design::baseline();
    let mut cases = Vec::new();
    for (n_fp, a_short) in [(3usize, 0.05), (3, 0.04), (6, 0.05)] {
        for (shift, shift_z) in [(0.0, 0.0), (0.01, -0.01), (0.02, -0.02), (-0.01, 0.01), (-0.02, 0.02), (0.0, 0.02), (0.0, -0.02)] {
            cases.push(
                RotatingEllipse {
                    n_fp,
                    coils_per_period: 18 / n_fp,
                    a_short,
                    shift,
                    shift_z,
                    ..base
                }
                .with_field(design::B0, design::R0),
            );
        }
    }
    let t0 = Instant::now();
    let rows: Vec<String> = cases
        .par_iter()
        .map(|c| {
            let f = design::filaments(c);
            let tr = design::tracer(c);
            let fps = tr.fixed_points(&f, (c.rc, 0.0));
            let n_x = fps.iter().filter(|p| !p.elliptic()).count();
            let head = format!(
                "n_fp {} x{:<2} a {:.2}/{:.2} shift R {:+.2} Z {:+.2}:",
                c.n_fp, c.coils_per_period, c.a_long, c.a_short, c.shift, c.shift_z
            );
            let Some(ax) = fps.iter().find(|p| p.elliptic()) else {
                return format!("{head} no elliptic axis ({} fixed points)", fps.len());
            };
            let axis = ax.rz;
            let path = tr.axis_path(&f, axis);
            let mut s = format!(
                "{head} axis R {:.4} ({} O, {n_x} X) iota0 {:.4} |",
                axis.0,
                fps.len() - n_x,
                ax.iota(c.n_fp)
            );
            for d in [0.005, 0.02, 0.035, 0.05] {
                let l = tr.trace(&f, (axis.0 + d, axis.1), &path, 60);
                match l.lost_after {
                    Some(p) => s += &format!("  {:>4.0}mm lost@{p}", d * 1e3),
                    None => s += &format!("  {:>4.0}mm {:.4}", d * 1e3, l.iota(c.n_fp).0),
                }
            }
            s
        })
        .collect();
    for r in rows {
        println!("{r}");
    }
    println!("{:.1} s", t0.elapsed().as_secs_f64());
}
