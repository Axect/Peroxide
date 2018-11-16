extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let a = poly(c!(1,3,2));
    let b = poly(c!(1,2));
    (a * b).print();
}