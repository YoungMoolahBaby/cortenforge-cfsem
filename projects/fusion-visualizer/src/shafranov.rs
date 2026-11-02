//! `shafranov`: the free-boundary solver against Shafranov's vertical-field formula.
//!
//! A circular plasma at major radius R0 needs an external vertical field
//! Bv = -(mu0 Ip / (4 pi R0)) [ln(8 R0 / a) + beta_p + li/2 - 3/2] (large aspect ratio). Here a
//! far coil pair supplies a nearly uniform vertical field, its current set each iteration so
//! the flux is equal at the inboard and outboard limiter points (R0 -+ a, 0): the plasma sits
//! centred at R0 with minor radius a. The field the pair then makes at R0 is compared with the
//! formula, using beta_p and li measured from the solution. The formula's error is O(a/R0), so
//! agreement should improve as a shrinks.

use std::f64::consts::PI;
use std::time::Instant;

use cfsem::MU_0;

use crate::equilibrium::{
    Coil, Control, Machine, Options, Profile, coil_filaments, filament_field, filament_flux, solve,
    summarize,
};
use crate::gs::{Grid, Solver};

const R0: f64 = 1.0;
const IP: f64 = 200e3;

fn vf_pair(current: f64) -> Vec<Coil> {
    let c = |z| Coil {
        name: "VF",
        r: 3.0,
        z,
        w: 0.1,
        h: 0.1,
        n: 3,
        current,
    };
    vec![c(-2.0), c(2.0)]
}

pub fn run() {
    let grid = Grid::new((0.4, 1.6), 121, (-0.6, 0.6), 121);
    let t0 = Instant::now();
    let solver = Solver::new(&grid, true);
    println!(
        "Grid {} x {} (dR {:.4} m), order-4 Delta* factorized in {:.2} s; Ip {:.0} kA, R0 {R0} m",
        grid.nr(),
        grid.nz(),
        grid.dr(),
        t0.elapsed().as_secs_f64(),
        IP / 1e3
    );
    let unit = coil_filaments(&vf_pair(1.0));
    let bz_per_amp = filament_field(&unit, &[R0], &[0.0]).1[0];
    let h = 1e-3;
    let bz = |r: f64| filament_field(&unit, &[r], &[0.0]).1[0];
    println!(
        "VF pair at R 3.0 m, Z +-2.0 m: {bz_per_amp:.3e} T/A at (R0, 0), decay index there {:.3}",
        -R0 / bz(R0) * (bz(R0 + h) - bz(R0 - h)) / (2.0 * h)
    );
    let verbose = std::env::var("FV_DEBUG").is_ok();
    println!(
        "  a/R0  VF current    Bv(R0)       Shafranov    diff     beta_p  li     a_lcfs  R_geo   iters  time"
    );
    for a in [0.3, 0.2, 0.1] {
        let t0 = Instant::now();
        let region: Vec<(f64, f64)> = (0..128)
            .map(|k| {
                let t = 2.0 * PI * k as f64 / 128.0;
                (R0 + 1.5 * a * t.cos(), 1.5 * a * t.sin())
            })
            .collect();
        let machine = Machine {
            coils: vf_pair(0.0),
            region,
            limiter: vec![(R0 - a, 0.0), (R0 + a, 0.0), (R0, 1.4 * a), (R0, -1.4 * a)],
            control: Some(Control {
                coils: vec![0, 1],
                points: [(R0 - a, 0.0), (R0 + a, 0.0)],
            }),
        };
        let prof = Profile {
            ip: IP,
            beta0: 0.3,
            a1: 1.0,
            a2: 2.0,
            r0: R0,
        };
        let opt = Options {
            verbose,
            ..Options::default()
        };
        let eq = solve(&grid, &solver, &machine, &prof, &opt, None);
        let s = summarize(&grid, &eq, &prof);
        let i_vf = eq.control_current.expect("controlled");
        let bv = bz_per_amp * i_vf;
        let formula = -MU_0 * s.ip / (4.0 * PI * s.r_geo)
            * ((8.0 * s.r_geo / s.a).ln() + s.beta_p + 0.5 * s.li - 1.5);
        println!(
            "  {:.2}  {:>10.1} A  {:+.4e}  {:+.4e}  {:+.2}%  {:.3}   {:.3}  {:.4}  {:.4}  {:>5}  {:.1} s{}",
            a / R0,
            i_vf,
            bv,
            formula,
            (bv / formula - 1.0) * 100.0,
            s.beta_p,
            s.li,
            s.a,
            s.r_geo,
            eq.iterations,
            t0.elapsed().as_secs_f64(),
            if eq.converged {
                ""
            } else {
                "  (not converged)"
            }
        );
        if (a - 0.2).abs() < 1e-9 {
            // The plasma flux on the grid against cfsem's Green's function for the same current,
            // at vacuum nodes away from the plasma.
            let (rn, zn): (Vec<f64>, Vec<f64>) = grid.nodes().map(|n| (n.0, n.1)).unzip();
            let da = grid.dr() * grid.dz();
            let src: Vec<(f64, f64, f64)> = (0..grid.len())
                .filter(|&k| eq.jphi[k] != 0.0)
                .map(|k| (rn[k], zn[k], eq.jphi[k] * da))
                .collect();
            let far: Vec<usize> = (0..grid.len())
                .filter(|&k| (rn[k] - R0).hypot(zn[k]) > a + 5.0 * grid.dr())
                .filter(|&k| !grid.is_boundary(k / grid.nz(), k % grid.nz()))
                .collect();
            let (rf, zf): (Vec<f64>, Vec<f64>) = far.iter().map(|&k| (rn[k], zn[k])).unzip();
            let green = filament_flux(&src, &rf, &zf);
            let (mut e, mut m) = (0.0_f64, 0.0_f64);
            for (n, &k) in far.iter().enumerate() {
                e = e.max((eq.psi_plasma[k] - green[n]).abs());
                m = m.max(green[n].abs());
            }
            println!(
                "    plasma flux, grid vs cfsem Green's sum at {} vacuum nodes: rel err {:.2e}; \
                 Ip from J {:.1} kA",
                far.len(),
                e / m,
                s.ip / 1e3
            );
        }
    }
}
