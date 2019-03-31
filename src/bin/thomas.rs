extern crate peroxide;
use peroxide::*;

fn main() {
    let eps: f64 = 1f64;
    let h: f64 = 0.001;
    let n: usize = 1000;
    let x = seq(0, 1, h);

    // Traditional
    //let a_init = numerical(1f64, n, (1f64, 3f64));
    //let a = cbind(x.to_matrix(), a_init);
    //a.print();
    //a.write_pickle("example_data/traditional.pickle").expect("Can't write pickle file");
    
    // Thomas
    let prev = vec![(eps + h/2f64)/h.powi(2); n-1];
    let center = vec![-2f64*eps/h.powi(2); n];
    let post = vec![(eps - h/2f64)/h.powi(2); n-1];
    let y = b(eps, h, n, (1f64, 3f64)).col(0);

    thomas(prev,center,post,y).print();
}

fn thomas(a: Vec<f64>, b: Vec<f64>, c: Vec<f64>, y: Vec<f64>) -> Matrix {
    let n = b.len();
    assert_eq!(a.len(), n-1);
    assert_eq!(c.len(), n-1);
    assert_eq!(y.len(), n);

    let mut gamma: Vec<f64> = Vec::new();
    let mut beta: Vec<f64> = Vec::new();

    let mut gam = 0f64;
    let mut bet = 0f64;

    for i in 0 .. n-1 {
        let p = a[i];
        let q = b[i];
        let r = c[i];
        let s = y[i];

        bet = (s - p*bet)/(p*gam + q);
        gam = -r/(p*gam + q);

        beta.push(bet);
        gamma.push(gam);
    }

    let mut result = vec![0f64; n];
    let mut x = 0f64;
    let mut idx = n-1;

    while gamma.len() > 0 {
        let gam_i = gamma.pop().unwrap();
        let bet_i = beta.pop().unwrap();

        x = gam_i * x + bet_i;
        result[idx] = x;
        idx -= 1;
    }
    result.to_matrix()
}


/// Differential equation solver
#[allow(non_snake_case)]
fn numerical(eps: f64, N: usize, bcs: (f64, f64)) -> Matrix {
    let h = 1f64 / (N as f64);
    let size = (N - 1, N - 1);

    let A = A(eps, h, size);
    let b = b(eps, h, N - 1, bcs);

    let x_inner = A.inv().unwrap() % b.clone();

    let mut x_data = cat(bcs.0, x_inner.data);
    x_data.push(bcs.1);

    matrix(x_data, N + 1, 1, x_inner.shape)
}

/// Differential operator
#[allow(non_snake_case)]
fn A(eps: f64, h: f64, size: (usize, usize)) -> Matrix {
    let center = -2f64 * eps / h.powi(2);
    let prev = (eps + h / 2f64) / h.powi(2);
    let post = (eps - h / 2f64) / h.powi(2);

    let row = size.0;
    let col = size.1;

    let mut m = zeros(row, col);

    for i in 0..row {
        for j in 0..col {
            let idx = (i, j);
            if i == (j + 1) {
                m[idx] = prev;
            } else if i == j {
                m[idx] = center;
            } else if j == (i + 1) {
                m[idx] = post;
            }
        }
    }
    m
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

