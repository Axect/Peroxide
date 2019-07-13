extern crate peroxide;
use peroxide::*;

fn main() {
    let a = ml_matrix("1 2;3 4");
    let pqlu = a.lu().unwrap();
    pqlu.l.print();
    pqlu.u.print();
    println!("{:?}", pqlu.p);
    println!("{:?}", pqlu.q);

    a.inv().unwrap().print();

    let b = ml_matrix("1 2 2;4 5 6;7 8 9");
    b.inv().unwrap().print();
}
