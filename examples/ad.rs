#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = AD::from_array([3f64, 1f64, 0f64]);
    println!("{:?}", a);
    println!("{:?}", &a * &a);
}
