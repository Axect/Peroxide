use crate::numerical::utils::jacobian;
use crate::structure::matrix::*;
use crate::traits::{
    mutable::MutFP,
    math::{Vector, Normed, Norm},
    num::{Number, NumberVector},
};

/// Newton-Raphson Method
pub fn newton<F>(init_cond: Vec<f64>, f: F, rtol: f64) -> Vec<f64>
where
    F: Fn(Vec<Number>) -> Vec<Number> + Copy,
{
    let mut x_next = init_cond;
    let mut x = x_next.clone();
    update(&mut x_next, f);
    let mut err = (x_next.sub_vec(&x)).norm(Norm::L2);

    while err >= rtol {
        x = x_next.clone();
        update(&mut x_next, f);
        err = (x_next.sub_vec(&x)).norm(Norm::L2);
    }

    x_next
}

fn update<F>(xs: &mut Vec<f64>, f: F)
where
    F: Fn(Vec<Number>) -> Vec<Number> + Copy,
{
    let j = jacobian(f, &xs);
    let pinv_j = j.pseudo_inv().unwrap();
    let fx = f(NumberVector::from_f64_vec(xs.clone())).to_f64_vec();

    xs.mut_zip_with(|x,y| x - y, &(pinv_j * fx))
}
