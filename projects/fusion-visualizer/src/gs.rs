//! Grad-Shafranov plumbing around cfsem's Delta* operator: a regular (R, Z) grid, the sparse
//! operator, and a factorized solve.
//!
//! cfsem works in total poloidal flux, Psi = 2 pi R A_phi (Wb), so the equation is
//! Delta* Psi = -2 pi mu0 R J_phi. Its operator matrix indexes node (i, j) as i * nz + j (R outer)
//! and has identity rows on the boundary, so boundary values go straight into the right side.

use cfsem::MU_0;
use cfsem::physics::gradshafranov::{gs_operator_order2, gs_operator_order4};
use faer::prelude::*;
use faer::sparse::{SparseColMat, Triplet};
use std::f64::consts::PI;

/// A regular (R, Z) grid.
#[derive(Clone, Debug)]
pub struct Grid {
    pub rs: Vec<f64>,
    pub zs: Vec<f64>,
}

impl Grid {
    pub fn new(r: (f64, f64), nr: usize, z: (f64, f64), nz: usize) -> Self {
        let lin = |a: f64, b: f64, n: usize| -> Vec<f64> {
            (0..n)
                .map(|k| a + (b - a) * k as f64 / (n - 1) as f64)
                .collect()
        };
        Grid {
            rs: lin(r.0, r.1, nr),
            zs: lin(z.0, z.1, nz),
        }
    }
    pub fn nr(&self) -> usize {
        self.rs.len()
    }
    pub fn nz(&self) -> usize {
        self.zs.len()
    }
    pub fn len(&self) -> usize {
        self.nr() * self.nz()
    }
    pub fn dr(&self) -> f64 {
        self.rs[1] - self.rs[0]
    }
    pub fn dz(&self) -> f64 {
        self.zs[1] - self.zs[0]
    }
    pub fn idx(&self, i: usize, j: usize) -> usize {
        i * self.nz() + j
    }
    pub fn is_boundary(&self, i: usize, j: usize) -> bool {
        i == 0 || j == 0 || i + 1 == self.nr() || j + 1 == self.nz()
    }
    /// (R, Z, flat index) for every node.
    pub fn nodes(&self) -> impl Iterator<Item = (f64, f64, usize)> + '_ {
        (0..self.nr()).flat_map(move |i| {
            (0..self.nz()).map(move |j| (self.rs[i], self.zs[j], self.idx(i, j)))
        })
    }
    /// A field sampled at every node.
    pub fn sample(&self, f: impl Fn(f64, f64) -> f64) -> Vec<f64> {
        self.nodes().map(|(r, z, _)| f(r, z)).collect()
    }
}

/// cfsem's Delta* operator as (vals, rows, cols) triplets.
pub fn operator(grid: &Grid, order4: bool) -> (Vec<f64>, Vec<usize>, Vec<usize>) {
    if order4 {
        gs_operator_order4(&grid.rs, &grid.zs)
    } else {
        gs_operator_order2(&grid.rs, &grid.zs)
    }
}

/// y = A x for triplets (duplicates summed).
pub fn apply(op: &(Vec<f64>, Vec<usize>, Vec<usize>), x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; x.len()];
    for ((v, r), c) in op.0.iter().zip(&op.1).zip(&op.2) {
        y[*r] += v * x[*c];
    }
    y
}

/// A factorized Delta* for repeated solves on one grid.
pub struct Solver {
    lu: faer::sparse::linalg::solvers::Lu<usize, f64>,
    n: usize,
}

impl Solver {
    pub fn new(grid: &Grid, order4: bool) -> Self {
        let (vals, rows, cols) = operator(grid, order4);
        let n = grid.len();
        let trips: Vec<Triplet<usize, usize, f64>> = vals
            .iter()
            .zip(&rows)
            .zip(&cols)
            .map(|((&val, &row), &col)| Triplet { row, col, val })
            .collect();
        let mat = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &trips)
            .expect("operator triplets are in range");
        let lu = mat.sp_lu().expect("Delta* factorizes");
        Solver { lu, n }
    }

    /// Psi for a right-hand side that already carries the boundary values.
    pub fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        let b = Col::<f64>::from_fn(self.n, |i| rhs[i]);
        let x = self.lu.solve(&b);
        (0..self.n).map(|i| x[i]).collect()
    }
}

/// Right-hand side for a toroidal current density (A/m^2) at the nodes, with Psi on the boundary
/// taken from `boundary(r, z)`.
pub fn rhs(grid: &Grid, jphi: &[f64], boundary: impl Fn(f64, f64) -> f64) -> Vec<f64> {
    (0..grid.nr())
        .flat_map(|i| (0..grid.nz()).map(move |j| (i, j)))
        .map(|(i, j)| {
            let (r, z) = (grid.rs[i], grid.zs[j]);
            if grid.is_boundary(i, j) {
                boundary(r, z)
            } else {
                -2.0 * PI * MU_0 * r * jphi[grid.idx(i, j)]
            }
        })
        .collect()
}

/// (Br, Bz) from Psi on the grid by 4th-order central differences (interior nodes at least two
/// from the edge; others are NaN): Br = -dPsi/dZ / (2 pi R), Bz = dPsi/dR / (2 pi R).
pub fn field_from_flux(grid: &Grid, psi: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let (nr, nz) = (grid.nr(), grid.nz());
    let (mut br, mut bz) = (vec![f64::NAN; grid.len()], vec![f64::NAN; grid.len()]);
    let d4 =
        |m2: f64, m1: f64, p1: f64, p2: f64, h: f64| (m2 - 8.0 * m1 + 8.0 * p1 - p2) / (12.0 * h);
    for i in 2..nr - 2 {
        for j in 2..nz - 2 {
            let at = |a: usize, b: usize| psi[grid.idx(a, b)];
            let dr = d4(
                at(i - 2, j),
                at(i - 1, j),
                at(i + 1, j),
                at(i + 2, j),
                grid.dr(),
            );
            let dz = d4(
                at(i, j - 2),
                at(i, j - 1),
                at(i, j + 1),
                at(i, j + 2),
                grid.dz(),
            );
            let k = 2.0 * PI * grid.rs[i];
            br[grid.idx(i, j)] = -dz / k;
            bz[grid.idx(i, j)] = dr / k;
        }
    }
    (br, bz)
}
