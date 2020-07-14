extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = AD1::new(2f64, 1f64);
    a.print();
    let a2 = AD2::from(a);

    let b = AD2::new(4f64, 4f64, 2f64);
    b.print();

    (a + b).print();
    (a - b).print();
    (a * b).print();
    (a / b).print();
    (a2 / b).print();

    b.ln().print();
    let c = AD2::default();
    c.iter().rev().collect::<AD2>().print();
}
