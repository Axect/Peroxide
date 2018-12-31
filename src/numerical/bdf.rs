/// Backward Differentiation Formula
/// Author: Axect

use structure::matrix::*;
use structure::vector::*;
use structure::dual::*;
use numerical::utils::*;
use util::non_macro::*;

/// BDF1 (Backward Euler Method)
///
/// xs_init: Non-autonomous data: `xs = (t, ys)`
/// h: Step size
/// rtol: Relative tolerance (1e-15 is recommended)
pub fn bdf1<F>(xs_init: Vec<f64>, f: F, h: f64, rtol: f64, num: usize) -> Matrix
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    let mut xs = xs_init.clone();
    let mut records = zeros(num + 1, xs_init.len());
    for i in 0 .. (num+1) {
        records.subs_row(i, xs.clone());
        xs = one_step_bdf1(xs.clone(), f, h, rtol);
    }

    records
}


/// One step for BDF1 (Backward Euler)
pub fn one_step_bdf1<F>(xs: Vec<f64>, f: F, h: f64, rtol: f64) -> Vec<f64>
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    // len(ys) = n where xs = (t, ys)
    let n = xs.len() - 1;

    // Decompose xs to (t, ys)
    let t = xs[0];
    let ys = xs.clone().into_iter()
        .skip(1)
        .collect::<Vec<f64>>();
    
    // new t = t + h
    let new_t = t + h;

    // One step forward euler for initial guess
    // y_{n+1}^{(0)} = y_n + hf(t_n, y_n)
    let new_ys_0 = ys.add(
        &(f(xs.conv_dual()).values().fmap(|y| h*y))
    );

    let mut iter_prev_ys = new_ys_0;
    let mut iter_next_ys = non_auto_update(new_t, ys.clone(), iter_prev_ys.clone(), f, h);
    let mut err = iter_next_ys.sub(&iter_prev_ys).norm();

    while err >= rtol {
        iter_prev_ys = iter_next_ys.clone();
        iter_next_ys = non_auto_update(new_t, ys.clone(), iter_prev_ys.clone(), f, h);
        err = iter_next_ys.sub(&iter_prev_ys).norm();
    }

    let new_xs = concat(vec![new_t], iter_next_ys);

    new_xs
}

#[allow(non_snake_case)]
fn non_auto_update<F>(t: f64, yn: Vec<f64>, ys: Vec<f64>, f: F, h: f64) -> Vec<f64>
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    let n = ys.len();
    let xs = concat(vec![t], ys.clone());
    let Df = non_auto_jacobian(xs.clone(), f);
    let I = eye(n);
    let DF = I - h * Df;
    let f_xs = f(xs.conv_dual()).values();

    // F(y_{n+1}^i) = y_{n+1}^i - y_n - hf(t_n+h, y_{n+1}^i) 
    let F_ys = ys.sub(&yn).sub(&(f_xs.fmap(|x| h*x)));

    ys.sub(&(DF.inv().unwrap() % F_ys).col(0))
}
