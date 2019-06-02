extern crate peroxide;
use peroxide::*;

fn main() {
    // t=1, x1=2, x2=3
    let x0 = c!(1, 2, 3);

    // jacobian at t=1, x1=2, x2=3
    jacobian(x0, f).print();
}

fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let t = xs[0];
    let x1 = xs[1];
    let x2 = xs[2];

    vec![x1 * (x2 * t).exp()]
}