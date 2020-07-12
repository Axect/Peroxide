#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = AD1::new(2f64, 1f64);
    a.print();

    let b = AD2::new(4f64, 4f64, 2f64);
    b.print();

    (a + b).print();
    (a - b).print();
}
