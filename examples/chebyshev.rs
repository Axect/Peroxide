extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let mut c = chebyshev_polynomial(0, SpecialKind::First);
    c.print();
    for i in 1 .. 11 {
        c = chebyshev_polynomial(i, SpecialKind::First);
        c.print();
    }
}
