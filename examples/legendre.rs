extern crate peroxide;
use peroxide::*;

fn main() {
    let a = gauss_legendre_quadrature(|x| 1f64 / x, 15, (1f64, 2f64));
    a.print();
}
