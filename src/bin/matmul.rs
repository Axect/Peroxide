extern crate peroxide;
use peroxide::*;

fn main() {
    let a = ml_matrix("1;2;3;4");
    a.print();
    let b = ml_matrix("1 2 3 4");
    b.print();

    (&b * &a).print();
    (&a.t() * &a).print();
}