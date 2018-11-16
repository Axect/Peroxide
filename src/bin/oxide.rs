extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let x = c!(-1,0,1);
    let y = c!(1,0,1);
    let l = lagrange_polynomial(x, y);
    l.print();
}