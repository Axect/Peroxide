extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let x = chebyshev_nodes(10, 0, 1);
    let y = x.fmap(|t| t.powf(3.) + 3. * t.powf(2.) - 1.);
    let l = lagrange_polynomial(x, y);
    l.print();

    let l2 = l.integral();
    l2.print();
    l2.eval(1).print();
    (l2.eval(1) - l2.eval(0)).print();

    let x2 = linspace!(0, 1, 10);
    let y2 = x2.fmap(|t| t.powf(3.) + 3. * t.powf(2.) - 1.);
    let r = lagrange_polynomial(x2, y2);
    r.print();

    let r2 = r.integral();
    r2.print();
    r2.eval(1).print();
    (r2.eval(1) - r2.eval(0)).print();
}