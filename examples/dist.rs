extern crate peroxide;
use peroxide::fuga::*;

const SAMPLE_SIZE: usize = 10000;

fn main() {
    example();

    println!("Hello, World!")
}

fn example() {
    let mu: f64 = 0.1;
    let b = Bernoulli(mu);
    let b_sample = b.sample(SAMPLE_SIZE);
    println!("Bernoulli at mu = {}", mu);
    println!("Number of samples: {}", SAMPLE_SIZE);
    println!("Theoretical mean: {}", b.mean());
    println!("Sample mean: {}", b_sample.mean());
    println!("Theoretical var: {}", b.var());
    println!("Sample var: {}", b_sample.var());
    println!("PDF at x = 1: {}", b.pdf(1));
    println!("PDF at x = 0: {}", b.pdf(0));
    println!("");

    let m = 0;
    let s = 1;
    let norm = Normal(m, s);
    let norm_sample = norm.sample(SAMPLE_SIZE);
    println!("Normal at (a,b) = ({},{})", m, s);
    println!("Number of samples: {}", SAMPLE_SIZE);
    println!("Theoretical mean: {}", norm.mean());
    println!("Sample mean: {}", norm_sample.mean());
    println!("Theoretical var: {}", norm.var());
    println!("Sample var: {}", norm_sample.var());
    println!("PDF at x = mean: {}", norm.pdf(0f64));
    println!("");

    let alpha = 3;
    let beta = 2;
    let be = Beta(alpha, beta);
    let be_sample = be.sample(SAMPLE_SIZE);
    println!("Beta at (a,b) = ({},{})", alpha, beta);
    println!("Number of samples: {}", SAMPLE_SIZE);
    println!("Theoretical mean: {}", be.mean());
    println!("Sample mean: {}", be_sample.mean());
    println!("Theoretical var: {}", be.var());
    println!("Sample var: {}", be_sample.var());
    println!("PDF at x = mean: {}", be.pdf(0.6));
    println!("PDF at x = mode: {}", be.pdf(2f64 / 3f64));
    println!("");

    let alpha2 = 4;
    let beta2 = 6;
    let g = Gamma(alpha2, beta2);
    let g_sample = g.sample(SAMPLE_SIZE);
    println!("Gamma at (a,b) = ({},{})", alpha2, beta2);
    println!("Number of samples: {}", SAMPLE_SIZE);
    println!("Theoretical mean: {}", g.mean());
    println!("Sample mean: {}", g_sample.mean());
    println!("Theoretical var: {}", g.var());
    println!("Sample var: {}", g_sample.var());
    println!("PDF at x = mean: {}", g.pdf(0.666666));
    println!("PDF at x = mode: {}", g.pdf(0.5));
    println!("");

    let m = 1.0;
    let s = 1.0;
    let lognorm = LogNormal(m, s);
    let lognorm_sample = lognorm.sample(SAMPLE_SIZE);
    println!("LogNormal at (m,s) = ({},{})", m, s);
    println!("Number of samples: {}", SAMPLE_SIZE);
    println!("Theoretical mean: {}", lognorm.mean());
    println!("Sample mean: {}", lognorm_sample.mean());
    println!("Theoretical var: {}", lognorm.var());
    println!("Sample var: {}", lognorm_sample.var());
    println!("PDF at x = mean: {}", lognorm.pdf(2.718281828459045));
    println!("PDF at x = mode: {}", lognorm.pdf(1.0));
    println!("");
}
