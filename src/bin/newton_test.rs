extern crate peroxide;
use peroxide::*;

fn main() {
    let init_vec = c!(1,2,3);
    newton(init_vec, f, 1e-15).print();
}

fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let x1 = xs[0];
    let x2 = xs[1];
    let x3 = xs[2];

    vec![
        x1.pow(2) - 2.*x1 + x2.pow(2) - x3 + 1.,
        x1*x2.pow(2) - x1 - 3.*x2 + x2*x3 + 2.,
        x1*x3.pow(2) - 3.*x3 + x2*x3.pow(2) + x1*x2,
    ]
}