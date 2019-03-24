extern crate peroxide;
use peroxide::*;

fn main() {
    let a = matrix(seq(1, 441, 1), 21, 21, Col);
    a.print();
}