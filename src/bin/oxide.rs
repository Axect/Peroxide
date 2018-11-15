extern crate peroxide;

use peroxide::*;
#[allow(unused_imports)]
use std::process;

#[allow(unused_must_use)]
fn main() {
    let a = linspace!(0,1,11);
    a.print();
    let b = seq!(0,1,0.1);
    b.print();
    let r1 = rand!();
    let r2 = rand!(4, 2);
    println!("{}", r1);
    r2.print();

    let v_u32 = Uniform::new(1u32, 11);
    v_u32.sample(10).print();
    let v_f64 = Uniform::new(1f64, 11f64);
    v_f64.sample(10).print();

    println!("{}", erf(1.0));

    let v_n = Normal::new(0, 1);
    v_n.sample(10).print();
    v_n.sample(1000).mean().print();
    v_n.sample(1000).sd().print();

    let c1 = c!(1,2,3);
    let c2 = c!(3,2,1);
    c1.add(&c2).print();
}