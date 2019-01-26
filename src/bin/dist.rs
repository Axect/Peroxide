extern crate peroxide;
use peroxide::*;

fn main() {
    let a = Normal(0, 1);
    a.sample(10).print();
    a.pdf(0).print();

    let b = Bernoulli(0.1);
    b.sample(10).print();
    b.pdf(0).print();
}