extern crate peroxide;
use peroxide::*;
use std::f64::consts::PI;

fn main() {
    // x^2 * x (at x = 1)
    let a = hyper_dual(1, 2, 2);
    let b = hyper_dual(1, 1, 0);
    (a * b).print();

    // x at PI
    let x = hyper_dual(PI, 1f64, 0f64);
    // sin(x)
    println!("Sin: ");
    x.sin().print();
    // cos(x)
    println!("Cos: ");
    x.cos().print();
    // tan(x)
    println!("Tan: ");
    x.tan().print();

    // e^(-x^2)
    // 1. y = x^2 at x=1
    let y = hyper_dual(1, 2, 2);
    // 2. e^(-y)
    (-y).exp().print();

    // ln(exp(x^2))
    y.exp().ln().print();

    // ln(x^2)
    y.ln().print();

    y.powi(2).print();
    y.sqrt().print();

    let z = hyper_dual(1, 1, 0);
    (z / y).print();
}
