extern crate peroxide;
use peroxide::*;

fn main() {
    #[cfg(feature = "specials")]
    let a = gauss_legendre_quadrature(|x| 1f64 / x, 15, (1f64, 2f64));
    #[cfg(feature = "specials")]
    a.print();
}
