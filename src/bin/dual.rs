extern crate peroxide;
use peroxide::*;

fn main() {
    let r = dual(3, 0);
    let m = dual(2, 1);
    (1f64 - m / r).print();
}
