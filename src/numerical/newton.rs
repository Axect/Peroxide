/// Newton-Raphson Method

use structure::matrix::*;
use structure::dual::*;
use structure::vector::*;
use util::print::*;
use numerical::utils::jacobian;

pub fn newton<F>(init_cond: Vec<f64>, f: F, rtol: f64) -> Vec<f64>
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    let mut x = init_cond;
    let mut x_next = update(x.clone(), f);
    let mut err = (x_next.sub(&x)).norm();

    while err >= rtol {
        x = x_next.clone();
        x_next = update(x_next.clone(), f);
        err = (x_next.sub(&x)).norm();
        x_next.print();
        err.print();
    }

    x_next
}

fn update<F>(xs: Vec<f64>, f: F) -> Vec<f64>
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    let j = jacobian(xs.clone(), f);
    let pinv_j = j.pseudo_inv().unwrap();
    let xs_dual = xs.conv_dual();
    let fx = f(xs_dual).values();

    xs.sub(&(pinv_j % fx).col(0))
}