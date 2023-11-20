extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = AD1(2f64, 1f64);
    a.print();
    let a2 = a.to_order(2);

    let b = AD2(4f64, 4f64, 2f64);
    b.print();

    (a + b).print();
    (a - b).print();
    (a * b).print();
    (a / b).print();
    (a2 / b).print();

    a.exp().print();
    a.ln().print();
    b.exp().print();
    b.ln().print();
    a.powi(2).print();
    b.sqrt().print();
    b.cos().print();
    b.sin().print();
    b.cosh().print();
    b.sinh().print();

    f(a, b).print();

    (b + 1f64).print();
    (b - 1f64).print();
    (b * 2f64).print();
    (b / 2f64).print();
    assert_eq!(-(1f64 - b), b - 1f64);
    assert_eq!(2f64 * b, b * 2f64);
    assert_eq!(1f64 / (2f64 / b), b / 2f64);
}

fn f(a: AD, b: AD) -> AD {
    a + b
}
