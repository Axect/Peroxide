extern crate peroxide;

use peroxide::*;
use std::f64::consts::PI;

#[allow(unused_must_use)]
fn main() {
    let x = chebyshev_nodes(10, 0., PI);
    x.print();
    let y = x.fmap(|t| t.sin());
    let l = lagrange_polynomial(x, y);
    l.print();

    let l2 = l.integral();
    l2.print();
    (l2.eval(PI) - l2.eval(0)).print();

    let x2 = linspace!(0., PI, 10);
    x2.print();
    let y2 = x2.fmap(|t| t.sin());
    let r = lagrange_polynomial(x2, y2);
    r.print();

    let r2 = r.integral();
    r2.print();
    (r2.eval(PI) - r2.eval(0)).print();
}