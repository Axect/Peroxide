extern crate peroxide;
use peroxide::*;

pub fn main() {
    let xs = c!(1, 1);
    jacobian(xs, f).print();

    let t = dual(1, 1);
    let x = dual(1, 0);

    f(vec![t, x]).print();
}

fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let t = xs[0];
    let x = xs[1];

    vec![t.pow(2) * 3. * x]
}
