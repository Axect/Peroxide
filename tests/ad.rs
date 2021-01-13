extern crate peroxide;
#[allow(unused_imports)]
use peroxide::prelude::*;

#[test]
fn add_test() {
    let a = AD2(4f64, 4f64, 2f64);
    let b = AD1(2f64, 1f64);
    assert_eq!(a + b, AD2(6f64, 5f64, 2f64));
    println!("a     : {:?}", a);
    println!("b     : {:?}", b);
    println!("Neg a : {:?}", -a);
    println!("a + b : {:?}", a + b);
    println!("a - b : {:?}", a - b);
    println!("a * b : {:?}", a * b);
    println!("a / b : {:?}", a / b);
    println!("e^a   : {:?}", a.exp());
    println!("ln(a) : {:?}", a.ln());
    println!("a^2   : {:?}", a.powi(2));
    println!("a^0.5 : {:?}", a.powf(0.5));
    println!("a.sin : {:?}", a.sin());
    println!("a.cos : {:?}", a.cos());
    println!("a.tan : {:?}", a.tan());
    println!("a.sinh: {:?}", a.sinh());
    println!("a.cosh: {:?}", a.cosh());
    println!("a.tanh: {:?}", a.tanh());
    println!("a ^ b : {:?}", a.pow(b));
    a.pow(b).print();
}
