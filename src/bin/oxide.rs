extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let x = matrix!(1;5;1, 5, 1, Col);
    println!("{}", x);
    let y = matrix(c!(3.7, 4.2, 4.9, 5.7, 6.0), 5, 1, Col);
    println!("{}", lm(&x, &y));
    println!("{}", lm!(y ~ x));
}
