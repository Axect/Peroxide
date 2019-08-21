extern crate peroxide;
use peroxide::*;

fn main() {
    let a = rand(1000, 1000);
    (&a * &a).print();
}
