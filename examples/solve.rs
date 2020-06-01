extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = ml_matrix("1 -1; 1 1");
    let b = ml_matrix("1; 1");
    let opt_x = solve(&a, &b);
    opt_x.unwrap().print();
}