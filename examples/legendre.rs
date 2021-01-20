extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = integrate(|x| x.sin(), (1f64, 2f64), GaussLegendre(15));
    let b = integrate(|x| x.sin(), (1f64, 2f64), GaussLegendre(29));
    a.print();
    b.print();
}
