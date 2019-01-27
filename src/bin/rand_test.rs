extern crate rand;
extern crate peroxide;

use rand::prelude::*;
use rand::distributions::{Gamma, Normal};
use peroxide::*;

fn main() {
    let mut rng1 = thread_rng();
    let mut rng2 = thread_rng();

    let n = Normal::new(0f64, 1f64);
    println!("{}", rng1.sample(n));
    println!("{}", rng2.sample(n));

    let g = Gamma::new(4f64, 1f64/6f64);
    let mut v = vec![0f64; 10000];

    for i in 0 .. 10000 {
        v[i] = rng1.sample(g);
    }

    v.mean().print();
    v.var().print();
}