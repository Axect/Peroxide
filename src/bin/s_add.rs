extern crate peroxide;
use peroxide::*;

fn main() {
    let a = vec![0f64; 1000_0000];
    a.s_add(1f64).print();
}