#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = AD::from_array([3f64, 1f64, 0f64]);
    a.print();
    a.powi(2).print();
    (&a * &a).print();
    a.sqrt().print();
    a.sin().print();
    a.cos().print();
    a.tan().print();
    a.exp().print();
    a.ln().print();
    a.log2().print();
    a.log10().print();
    a.sinh().print();
    a.cosh().print();
    a.tanh().print();
}
