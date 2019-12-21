extern crate peroxide;
use peroxide::*;
use Number::D;

fn main() {
    let x = vec![D(dual(1f64, 1f64)), D(dual(2, 1))];
    jacobian(f, x.to_f64_vec()).print();
}

fn f(v: Vec<Number>) -> Vec<Number> {
    map(g, &v)
}

fn g(x: Number) -> Number {
    let d = x.to_dual();
    if d.value() < 2f64 {
        x.powi(2)
    } else {
        x.powi(3)
    }
}
