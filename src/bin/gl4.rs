extern crate peroxide;
use peroxide::*;

fn main() {
    let (k1, k2) = k_newton(f, 0., c!(1, 2), 1., 1e-15);
    k1.print();
    k2.print();
}

fn f(t: Dual, xs: Vec<Dual>) -> Vec<Dual> {
    let x = xs[0];
    let y = xs[1];

    vec![x + y, x - y]
}