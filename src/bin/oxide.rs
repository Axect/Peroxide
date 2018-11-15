extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let a = Polynomial::new(c!(1,3,2,5,4));
    a.print();
    let b = Polynomial::new(c!(3,2,1));
    (a.clone() - b).print();
    (a * 2).print();
}