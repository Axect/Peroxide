extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    // Declare polynomial
    let a = poly(c!(1,3,2));
    a.print();      // x^2 + 3x + 2
    a.eval(1);   // Evaluate when x = 1 -> 6.0

    let b = poly(c!(1,2,3,4)); // x^3 + 2x^2 + 3x + 4
    (a.clone() + b.clone()).print(); // x^3 + 3x^2 + 6x + 6
    (a.clone() - b.clone()).print(); // -x^3 - x^2 - 2
    (a.clone() * b.clone()).print();
}