extern crate peroxide;
use peroxide::*;

fn main() {
    let a = ml_matrix("1 2 5;4 5 6;7 8 9");
    a.det().print();
}