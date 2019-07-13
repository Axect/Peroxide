use numerical::utils::jacobian;
use structure::dual::*;
use structure::matrix::*;
use structure::vector::*;
use util::non_macro::{cat, concat, eye, zeros};

/// Value of 3f64.sqrt()
const SQRT3: f64 = 1.7320508075688772;

/// Butcher tableau for Gauss_legendre 4th order
const GL4: [[f64; 3]; 2] = [
    [0.5 - SQRT3 / 6f64, 0.25, 0.25 - SQRT3 / 6f64],
    [0.5 + SQRT3 / 6f64, 0.25 + SQRT3 / 6f64, 0.25],
];

/// Gauss-Legendre 4th order method
///
/// # Description
///
/// ## 1. Find `k1, k2`
/// * `k1 = f(t + p0*h, y + h(p1*k1 + p2*k2))`
/// * `k2 = f(t + q0*h, y + h(q1*k1 + q2*k2))`
///
/// ## 2. Iteration
/// * `y_{n+1} = y_n + 0.5*h*(k1 + k2)`
pub fn gl4<F>(f: F, t_init: f64, y_init: Vec<f64>, h: f64, rtol: f64, num: usize) -> Matrix
where
    F: Fn(Dual, Vec<Dual>) -> Vec<Dual> + Copy,
{
    let mut t = t_init;
    let mut y_curr = y_init.clone();
    let mut records = zeros(num + 1, y_curr.len() + 1);
    records.subs_row(0, cat(t, y_curr.clone()));

    for i in 0..num {
        let (k1, k2) = k_newton(f, t, y_curr.clone(), h, rtol);
        y_curr = y_curr.add(&k1.fmap(|x| 0.5 * x * h).add(&k2.fmap(|x| 0.5 * x * h)));
        t += h;
        records.subs_row(i + 1, cat(t, y_curr.clone()))
    }

    records
}

/// Newton's Method for find k in GL4
///
/// # Description
///
/// ## 0. Initial Guess by Euler method
/// * `k1 = f(t, y)`
/// * `k2 = f(t, y)`
///
/// ## 1. Combine below two equations to one equation
/// * `k1 = f(t1, y + h(p1*k1 + p2*k2))`
/// * `k2 = f(t2, y + h(q1*k1 + q2*k2))`
/// * `k = g(k)`
///
/// ## 2. Obtain Jacobian
/// * `DG(k^l) = I - Dg(k^l)`
///
/// ## 3. Iteration by Newton's Method
/// * `k^{l+1} = k^l - DG^{-1}G(k^l)`
#[allow(non_snake_case)]
pub fn k_newton<F>(f: F, t: f64, y: Vec<f64>, h: f64, rtol: f64) -> (Vec<f64>, Vec<f64>)
where
    F: Fn(Dual, Vec<Dual>) -> Vec<Dual> + Copy,
{
    let t1 = dual(t + GL4[0][0] * h, 0.);
    let t2 = dual(t + GL4[1][0] * h, 0.);
    let tn = dual(t, 0.);
    let yn = y.conv_dual();
    let n = y.len();

    // 0. Initial Guess
    let k1_init = f(tn, yn.clone()).values();
    let k2_init = f(tn, yn.clone()).values();
    let mut k_curr = concat(k1_init.clone(), k2_init.clone());
    let mut err = 1f64;

    // 1. Combine two functions to one function
    let g = |k: Vec<Dual>| -> Vec<Dual> {
        let k1 = k.take(n);
        let k2 = k.skip(n);
        concat(
            f(
                t1,
                yn.add(
                    &k1.fmap(|x| x * GL4[0][1] * h)
                        .add(&k2.fmap(|x| x * GL4[0][2] * h)),
                ),
            ),
            f(
                t2,
                yn.add(
                    &k1.fmap(|x| x * GL4[1][1] * h)
                        .add(&k2.fmap(|x| x * GL4[1][2] * h)),
                ),
            ),
        )
    };

    // 2. Obtain Jacobian
    let I = eye(2 * n);

    let mut Dg = jacobian(k_curr.clone(), g.clone());
    let mut DG = I.clone() - Dg.clone();
    let mut DG_inv = DG.inv().unwrap();
    let mut G = k_curr.sub(&g.clone()(k_curr.conv_dual()).values());
    let mut num_iter: usize = 0;

    // 3. Iteration
    while err >= rtol && num_iter <= 10 {
        let k_prev = k_curr.clone();
        let DGG = DG_inv.clone() * G.clone();
        k_curr = k_curr.sub(&DGG.col(0));
        Dg = jacobian(k_curr.clone(), g.clone());
        DG = I.clone() - Dg.clone();
        DG_inv = DG.inv().unwrap();
        G = k_curr.sub(&g.clone()(k_curr.conv_dual()).values());
        err = k_curr.sub(&k_prev).norm();
        num_iter += 1;
    }

    (k_curr.take(n), k_curr.skip(n))
}
