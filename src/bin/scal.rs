extern crate peroxide;
use peroxide::*;

fn main() {
    let a = vec![0f64; 1000_0000];
    let c = 1f64;
    a.s_mul(c).print();
}