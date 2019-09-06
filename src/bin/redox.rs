extern crate peroxide;
use peroxide::*;

fn main() {
    let a = RedoxVector::from_vec(vec![0f64; 1000_0000]);
    let b = RedoxVector::from_vec(vec![0f64; 1000_0000]);

    (*(a + b)).print();
}