extern crate peroxide;
use peroxide::*;

fn main() {
    let a = ml_matrix("1 2;3 4");
    let b = ml_matrix("5 6;7 8");

    let c = &a + &b;

    c.print();
    a.print();
    b.print();
}
