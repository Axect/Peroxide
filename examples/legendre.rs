extern crate peroxide;
use peroxide::*;

fn main() {
    let a = integrate(|x| 1f64 / x, (1f64, 2f64), GaussLegendre(15));
    a.print();
}
