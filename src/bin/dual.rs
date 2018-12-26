extern crate peroxide;
use peroxide::*;

fn main() {
    let x = c!(1,2,3);
    let dx = c!(1,1,1);
    let xs = merge_dual(x,dx);
    xs.print();
    xs.values().print();
    xs.slopes().print();

    let j = jacobian(c!(1, 1), f);
    j.print();
}

fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let x = xs[0];
    let y = xs[1];

    vec![
        x - y,
        x + 2.*y,
    ]
}