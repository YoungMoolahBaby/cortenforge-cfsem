//! `shape`: a small tokamak whose plasma shape follows its PF coil currents.
//!
//! R0 = 1 m, Ip = 200 kA, an elliptical vessel. The vertical-field pair holds the radial position
//! by feedback (equal flux at R0 -+ 0.35 m on the midplane); the elongation pair (above and
//! below the vessel) and the triangularity pair (by the centre column) are set by hand, and the
//! shape is read off the last closed flux surface. Up-down symmetric cases only: the vertical
//! position is not controlled, and elongated plasmas are vertically unstable.

use std::f64::consts::PI;
use std::time::Instant;

use crate::equilibrium::{Coil, Control, Equilibrium, Machine, Options, Profile, solve, summarize};
use crate::gs::{Grid, Solver};
use crate::plot::equilibrium_svg;

pub const R0: f64 = 1.0;
pub const IP: f64 = 200e3;

/// The machine with elongation and triangularity pair currents (A-turns per coil).
pub fn machine(el: f64, tr: f64) -> Machine {
    let mut coils = vec![
        Coil {
            name: "VF",
            r: 2.0,
            z: 1.2,
            w: 0.15,
            h: 0.15,
            n: 3,
            current: 0.0,
        },
        Coil {
            name: "VF",
            r: 2.0,
            z: -1.2,
            w: 0.15,
            h: 0.15,
            n: 3,
            current: 0.0,
        },
        Coil {
            name: "EL",
            r: 1.0,
            z: 1.05,
            w: 0.3,
            h: 0.1,
            n: 4,
            current: el,
        },
        Coil {
            name: "EL",
            r: 1.0,
            z: -1.05,
            w: 0.3,
            h: 0.1,
            n: 4,
            current: el,
        },
        Coil {
            name: "TR",
            r: 0.3,
            z: 0.55,
            w: 0.1,
            h: 0.1,
            n: 3,
            current: tr,
        },
        Coil {
            name: "TR",
            r: 0.3,
            z: -0.55,
            w: 0.1,
            h: 0.1,
            n: 3,
            current: tr,
        },
    ];
    coils.push(Coil {
        name: "CS",
        r: 0.2,
        z: 0.0,
        w: 0.1,
        h: 1.4,
        n: 7,
        current: 0.0,
    });
    let vessel: Vec<(f64, f64)> = (0..160)
        .map(|k| {
            let t = 2.0 * PI * k as f64 / 160.0;
            (R0 + 0.42 * t.cos(), 0.68 * t.sin())
        })
        .collect();
    Machine {
        coils,
        limiter: vessel.clone(),
        region: vessel,
        control: Some(Control {
            coils: vec![0, 1],
            points: [(R0 - 0.35, 0.0), (R0 + 0.35, 0.0)],
        }),
    }
}

pub fn grid() -> Grid {
    Grid::new((0.4, 1.7), 131, (-0.9, 0.9), 181)
}

pub fn profile() -> Profile {
    Profile {
        ip: IP,
        beta0: 0.3,
        a1: 1.0,
        a2: 2.0,
        r0: R0,
    }
}

pub const CASES: [(&str, f64, f64); 6] = [
    ("baseline", 0.0, 0.0),
    ("elongation 40 kA", 40e3, 0.0),
    ("elongation 80 kA", 80e3, 0.0),
    ("elongation 120 kA", 120e3, 0.0),
    ("elongation 80 kA, triangularity +60 kA", 80e3, 60e3),
    ("elongation 80 kA, triangularity -60 kA", 80e3, -60e3),
];

pub fn run() {
    let grid = grid();
    let t0 = Instant::now();
    let solver = Solver::new(&grid, true);
    println!(
        "Small tokamak: R0 {R0} m, Ip {:.0} kA; grid {} x {}, order-4 Delta* factorized in {:.1} s",
        IP / 1e3,
        grid.nr(),
        grid.nz(),
        t0.elapsed().as_secs_f64()
    );
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("out");
    std::fs::create_dir_all(&out).expect("create out/");
    println!(
        "  case                                      VF (fb)    R_geo   a      kappa  d_up    d_low   beta_p  li     edge             iters"
    );
    let prof = profile();
    let mut prev: Option<Equilibrium> = None;
    for (k, (name, el, tr)) in CASES.iter().enumerate() {
        let m = machine(*el, *tr);
        let eq = solve(
            &grid,
            &solver,
            &m,
            &prof,
            &Options::default(),
            prev.as_ref(),
        );
        let s = summarize(&grid, &eq, &prof);
        println!(
            "  {name:<40}  {:>7.1} kA  {:.3}  {:.3}  {:.3}  {:+.3}  {:+.3}  {:.3}   {:.3}  {:<16} {:>4}{}",
            eq.control_current.unwrap_or(0.0) / 1e3,
            s.r_geo,
            s.a,
            s.kappa,
            s.delta_upper,
            s.delta_lower,
            s.beta_p,
            s.li,
            match eq.boundary {
                crate::equilibrium::Boundary::Limiter(r, z) => format!("wall {r:.2},{z:+.2}"),
                crate::equilibrium::Boundary::XPoint(r, z) => format!("X {r:.2},{z:+.2}"),
            },
            eq.iterations,
            if eq.converged { "" } else { " (not converged)" }
        );
        let mut placed = m;
        if let Some(i) = eq.control_current {
            placed.coils[0].current = i;
            placed.coils[1].current = i;
        }
        let path = out.join(format!("shape_{k}.svg"));
        equilibrium_svg(&grid, &eq, &s, &placed, name, &path);
        prev = Some(eq);
    }
    println!("  plots: {}", out.join("shape_*.svg").display());
}
