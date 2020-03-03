//! Lanczos approximation Coefficient generator

use std::f64::consts::PI;
use structure::matrix::Matrix;
use util::non_macro::zeros;
use util::useful::sgn;
use statistics::ops::{C, factorial, double_factorial};

const G: f64 = 5f64;
const N: usize = 7;

// Lanczos g=5, n=7
const LG5N7: [f64; 7] = [
    1.000000000189712,
    76.18009172948503,
    -86.50532032927205,
    24.01409824118972,
    -1.2317395783752254,
    0.0012086577526594748,
    -0.00000539702438713199
];

pub fn ln_gamma_approx(z: f64) -> f64 {
    let z = z - 1f64;
    let base = z + G + 0.5;
    let mut s = 0f64;
    for i in 1 .. N {
        s += LG5N7[i] / (z + i as f64);
    }
    s += LG5N7[0];
    (2f64 * PI).sqrt().ln() + s.ln() - base + base.ln() * (z + 0.5)
}

pub fn gamma_approx(z: f64) -> f64 {
    if z > 1f64 {
        let z_int = z as usize;
        if z - (z_int as f64) == 0f64 {
            return factorial(z_int-1) as f64;
        }
    }

    if z < 0.5 {
        PI / ((PI * z).sin() * gamma_approx(1f64 - z))
    } else {
        ln_gamma_approx(z).exp()
    }
}

/// Lanczos Approximation Coefficient
pub fn tlg1(g: f64, n: usize) -> Matrix {
    lanczos_coeff(g, n-1) * g.exp() / (2f64 * std::f64::consts::PI).sqrt()
}

fn lanczos_coeff(g: f64, n: usize) -> Matrix {
    let m = dr_gen(n) * b_gen(n) * (c_gen(n) * dc_gen(n));
    let f = f_gen(g, n);
    m * f
}

fn b_gen(n: usize) -> Matrix {
    Matrix::from_index(
        |i, j| {
            if i == 0 {
                1f64
            } else if j >= i {
                sgn(j - i) * C(i+j-1, j-i) as f64
            } else {
                0f64
            }
        },
        (n+1, n+1)
    )
}

fn c_gen(n: usize) -> Matrix {
    Matrix::from_index(
        |i, j| {
            if i == 0 && j == 0 {
                0.5
            } else if j > i {
                0f64
            } else {
                sgn(i-j) * 4f64.powi(j as i32) * (i as f64) * (C(i+j, 2*j) as f64) / (i+j) as f64
            }
        },
        (n+1, n+1)
    )
}

fn dc_gen(n: usize) -> Matrix {
    let mut m = zeros(n+1, n+1);
    m[(0,0)] = 2f64;
    for i in 1 .. n+1 {
        m[(i,i)] = 2f64 * double_factorial(2*i-1) as f64;
    }
    m
}

fn dr_gen(n: usize) -> Matrix {
    let mut m = zeros(n+1, n+1);
    m[(0,0)] = 1f64;
    for i in 1 .. n+1 {
        m[(i,i)] = - ((i * C(2*i-1, i)) as f64);
    }
    m
}

fn f(g: f64, n: usize) -> f64 {
    2f64.sqrt() * (n as f64 + 0.5).exp() / (2f64 * (n as f64 + g) + 1f64).powf(n as f64 + 0.5)
}

fn f_gen(g: f64, n: usize) -> Vec<f64> {
    let mut v = vec![0f64; n+1];
    for i in 0 .. n+1 {
        v[i] = f(g, i);
    }
    v
}