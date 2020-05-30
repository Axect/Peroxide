extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = rand(1000, 1000);
    (&a * &a).print();
}
