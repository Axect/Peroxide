extern crate peroxide;
use peroxide::*;

fn main() {
    let eps: f64 = 1f64;
    let h: f64 = 0.001;
    let n: usize = 1000;

    // Thomas
    let prev = vec![(eps + h / 2f64) / h.powi(2); n - 1];
    let center = vec![-2f64 * eps / h.powi(2); n];
    let post = vec![(eps - h / 2f64) / h.powi(2); n - 1];
    let y = b(eps, h, n, (1f64, 3f64)).col(0);

    let a = tdma(prev, center, post, y);
    a.print();
    a.write_single_pickle("example_data/tdma.pickle")
        .expect("Can't write pickle file");
}

/// RHS with Dirichlet Boundary Condition
fn b(eps: f64, h: f64, size: usize, bcs: (f64, f64)) -> Matrix {
    let first = -1f64 - 1f64 / h.powi(2) * (eps + h / 2f64) * bcs.0;
    let last = -1f64 - 1f64 / h.powi(2) * (eps - h / 2f64) * bcs.1;

    let mut col_mat = matrix(vec![-1f64; size], size, 1, Col);
    col_mat[(0, 0)] = first;
    col_mat[(size - 1, 0)] = last;
    col_mat
}
