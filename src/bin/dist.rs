extern crate peroxide;
use peroxide::*;

const sample_size: usize = 10000;

fn main() {
//    let a = Normal(0, 1);
//    let a_sample = a.sample(sample_size);
//    a_sample.mean().print();
//    a_sample.var().print();
//    a.mean().print();
//    a.var().print();
//    println!("");

    let mu: f64 = 0.1;
    let b = Bernoulli(mu);
    let b_sample = b.sample(sample_size);
    println!("Bernoulli at mu = {}",    mu);
    println!("Number of samples: {}",   sample_size);
    println!("Theoretical mean: {}",    b.mean());
    println!("Sample mean: {}",         b_sample.mean());
    println!("Theoretical var: {}",     b.var());
    println!("Sample var: {}",          b_sample.var());
    println!("PDF at x = 1: {}",        b.pdf(1));
    println!("PDF at x = 0: {}",        b.pdf(0));
    println!("");

    let alpha = 3;
    let beta = 2;
    let be = Beta(alpha, beta);
    let be_sample = be.sample(sample_size);
    println!("Beta at (a,b) = ({}, {})",    alpha, beta);
    println!("Number of samples: {}",       sample_size);
    println!("Theoretical mean: {}",        be.mean());
    println!("Sample mean: {}",             be_sample.mean());
    println!("Theoretical var: {}",         be.var());
    println!("Sample var: {}",              be_sample.var());
    println!("PDF at x = mean: {}",         be.pdf(0.6));
    println!("PDF at x = mode: {}",         be.pdf(2f64/3f64));
    println!("");
}