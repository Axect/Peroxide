use crate::numerical::utils::jacobian;
use crate::structure::matrix::*;
use crate::structure::ad::*;
use crate::traits::{
    math::{Norm, Normed, Vector},
    mutable::MutFP,
};

/// Newton-Raphson Method
pub fn newton<F: Fn(&Vec<AD>) -> Vec<AD> + Copy>(init_cond: Vec<f64>, f: F, rtol: f64) -> Vec<f64>
{
    let mut x_next = init_cond;
    let mut x = x_next.clone();
    update(&mut x_next, f);
    let mut err = x_next.sub_vec(&x).norm(Norm::L2);

    while err >= rtol {
        x = x_next.clone();
        update(&mut x_next, f);
        err = x_next.sub_vec(&x).norm(Norm::L2);
    }

    x_next
}

fn update<F: Fn(&Vec<AD>) -> Vec<AD> + Copy>(xs: &mut Vec<f64>, f: F)
{
    let j = jacobian(f, &xs);
    let pinv_j = j.pseudo_inv();
    //let fx = f(NumberVector::from_f64_vec(xs.clone())).to_f64_vec();
    let xs_ad: Vec<AD> = xs.iter().map(|&t| AD::from(t)).collect();
    let fx: Vec<f64> = f(&xs_ad).iter().map(|&t| t.into()).collect();

    xs.mut_zip_with(|x, y| x - y, &(pinv_j * fx))
}
