#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = ml_matrix("1 -1; 1 1");
    let b = c!(1, 1);
    let opt_x = a.solve(&b, LU);
    opt_x.print();
}