extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let r = dual(3, 0);
    let m = dual(2, 1);
    (1f64 - m / r).print();

    m.pow(r).print(); // (x^3, 3x^2) = (8, 12)
    r.pow(m).print(); // (3^x, 3^x * ln(3)) = (9, 9.8875...)
    m.pow(m).print(); // (x^x, x^x(ln(x) + 1)) = (4, 6.7725...)
}
