extern crate peroxide;

use peroxide::*;

#[allow(unused_must_use)]
fn main() {
    // Cubic spline example
    let x = c!(0.9, 1.3, 1.9, 2.1);
    let y = c!(1.3, 1.5, 1.85, 2.1);

    let s = cubic_spline(x, y);

    for i in 0 .. s.len() {
        s[i].print();
    }

    // Least square example
    let a = c!(1,2,3,4,5);
    let b = c!(1.2, 1.8, 3.2, 3.8, 5.0);
    let ls = least_square(a, b);
    ls.print();
    println!("{:?}", ls.coef);
}