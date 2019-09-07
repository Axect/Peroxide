extern crate peroxide;
use peroxide::*;

fn main() {
    let a = vec![0f64; 1000_0000];
    let b = vec![0f64; 1000_0000];

    a.dot(&b).print();
}