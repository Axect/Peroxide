extern crate peroxide;
use peroxide::*;

fn main() {
    let a = Normal(0, 1);
    a.sample(10).print();
    a.pdf(0).print();
}