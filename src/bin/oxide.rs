extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let a = c!(1,2,3,4,5).to_matrix();
    let b = a.clone() + Normal::new(0,1).sample(5).to_matrix();
    lm!(b ~ a).print();
}