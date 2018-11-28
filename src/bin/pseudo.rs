extern crate peroxide;
use peroxide::*;

fn main() {
    let a = matrix!(1;4;1, 2, 2, Row);
    a.pseudo_inv().unwrap().print();
    a.inv().unwrap().print();
}