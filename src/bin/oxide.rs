extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let a = poly(c!(3,2,1));
    a.diff().print();
    a.integral().print();
}