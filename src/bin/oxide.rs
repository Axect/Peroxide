extern crate peroxide;

use peroxide::*;

fn main() {
    let a = matrix!(1;4;1, 2, 2, Row);
    let pqlu = a.lu().unwrap();
    println!("{:?}\n{:?}\n{}\n{}", pqlu.p, pqlu.q, pqlu.l, pqlu.u);
    println!("{}", a.det());
    println!("{}", a.inv().unwrap());
}