extern crate peroxide;
use peroxide::*;

const sample_size: usize = 10000;

fn main() {
    let a = Normal(0, 1);
    let a_sample = a.sample(sample_size);
    a_sample.mean().print();
    a_sample.var().print();
    a.mean().print();
    a.var().print();
    println!("");

    let b = Bernoulli(0.1);
    b.sample(sample_size).mean().print();
    b.pdf(1).print();
    println!("");

    let be = Beta(3, 2);
    let be_sample = be.sample(sample_size);
    be_sample.mean().print();
    be_sample.var().print();
    be.mean().print();
    be.var().print();
}