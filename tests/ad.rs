extern crate peroxide;
#[allow(unused_imports)]
use peroxide::prelude::*;

#[test]
fn ad_test() {
    let a = AD2(4f64, 4f64, 2f64);
    let b = AD1(2f64, 1f64);
    let c = AD2(0.25f64, 1f64, 2f64);
    assert_eq!(a + b, AD2(6f64, 5f64, 2f64));
    println!("a      : {:?}", a);
    println!("b      : {:?}", b);
    println!("Neg a  : {:?}", -a);
    println!("a + b  : {:?}", a + b);
    println!("a - b  : {:?}", a - b);
    println!("a * b  : {:?}", a * b);
    println!("a / b  : {:?}", a / b);
    println!("e^a    : {:?}", a.exp());
    println!("ln(a)  : {:?}", a.ln());
    println!("a^2    : {:?}", a.powi(2));
    println!("a^0.5  : {:?}", a.powf(0.5));
    println!("a^b    : {:?}", a.pow(b));
    println!("a.sin  : {:?}", a.sin());
    println!("a.cos  : {:?}", a.cos());
    println!("a.tan  : {:?}", a.tan());
    println!("a.sinh : {:?}", a.sinh());
    println!("a.cosh : {:?}", a.cosh());
    println!("a.tanh : {:?}", a.tanh());
    println!("a ^ b  : {:?}", a.pow(b));
    println!("c.asin : {:?}", c.asin());
    println!("c.acos : {:?}", c.acos());
    println!("c.atan : {:?}", c.atan());
    println!("a.asinh: {:?}", a.asinh());
    println!("a.acosh: {:?}", a.acosh());
    println!("c.atanh: {:?}", c.atanh());
}
