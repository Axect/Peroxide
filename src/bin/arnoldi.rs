extern crate peroxide;
use peroxide::*;
use std::ops::Mul;
use std::cmp::max;

fn main() {
    let m = ml_matrix("1 2;3 4");
    let b = c![1, 0.5];
    let a = arnoldi_iteration(m, b, 5);
    a.print();
}


#[allow(dead_code)]
fn givens(alpha: f64, beta: f64) -> (f64, f64) {
    let mut tau = 0f64;
    let mut gamma = 0f64;
    let mut sigma = 0f64;

    if alpha == 0f64 || beta == 0f64 {
        tau = 0f64;
        let denominator = (alpha.powi(2) + beta.powi(2)).sqrt();
        gamma = alpha / denominator;
        sigma = beta / denominator;
    } else if alpha.abs() >= beta.abs() {
        tau = beta / alpha;
        gamma = 1f64 / (1f64 + tau.powi(2)).sqrt();
        sigma = gamma * tau;
    } else {
        tau = alpha / beta;
        sigma = 1f64 / (1f64 + tau.powi(2)).sqrt();
        gamma = sigma * tau;
    }

    (gamma, sigma)
}

#[allow(dead_code)]
fn arnoldi_iteration(a: Matrix, b: Vec<f64>, n: usize) -> Matrix {
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut x: Vec<Vec<f64>> = Vec::new();
    let mut y: Vec<Vec<f64>> = Vec::new();

    let q_temp = b.fmap(|t| t / b.norm());
    q.push(q_temp);

    // Make basis of Krylov space
    for k in 1 .. n {
        let x_temp = (a.clone() * q[k-1].clone()).data;
        let mut y_temp = x_temp.clone();
        for j in 0 .. k {
            let inner = x_temp.dot(&q[j]);
            y_temp = y_temp.sub(&q[j].fmap(|t| t * inner));
        }
        q.push(y_temp.fmap(|t| t / y_temp.norm()));
        x.push(x_temp); // n-1
        y.push(y_temp); // n-1
    }
    x.push((a.clone() * q[n-1].clone()).data);

    // Construct Upper Hessenberg matrix
    // h: n x (n-1)
    let mut h = zeros(n, n);
    let q0 = q[0].clone();
    for k in 0 .. n {
        h[(0, k)] = q0.dot(&x[k]);
    }
    for j in 1 .. n {
        let qj = q[j].clone();
        h[(j, j-1)] = y[j-1].norm();
        for k in j .. n {
            h[(j, k)] = qj.dot(&x[k]);
        }
    }
    h
}

#[allow(dead_code)]
fn normalized_vec(n: usize) -> Matrix {
    let v = rand(n, 1);
    v.fmap(|x| x / v.norm(Frobenius))
}
