use numerical::utils::jacobian;
use operation::number::{Number, NumberVector};
/// Newton-Raphson Method
use structure::matrix::*;
use structure::vector::*;

pub fn newton<F>(init_cond: Vec<f64>, f: F, rtol: f64) -> Vec<f64>
where
    F: Fn(Vec<Number>) -> Vec<Number> + Copy,
{
    let mut x = init_cond;
    let mut x_next = update(&x, f);
    let mut err = (x_next.sub(&x)).norm();

    while err >= rtol {
        x = x_next.clone();
        x_next = update(&x_next, f);
        err = (x_next.sub(&x)).norm();
    }

    x_next
}

fn update<F>(xs: &Vec<f64>, f: F) -> Vec<f64>
where
    F: Fn(Vec<Number>) -> Vec<Number> + Copy,
{
    let j = jacobian(f, &xs);
    let pinv_j = j.pseudo_inv().unwrap();
    let fx = f(NumberVector::from_f64_vec(xs.clone())).to_f64_vec();

    xs.sub(&(pinv_j * fx).col(0))
}
