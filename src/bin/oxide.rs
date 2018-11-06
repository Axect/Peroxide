extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let a = matrix!(1;4;1, 2, 2, Row);
    a.col(0).print();
    a.col(1).print();
    a.row(0).print();
    a.row(1).print();
}