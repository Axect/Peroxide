extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let x_f64 = 2f64;
    let x_ad1 = AD1(2f64, 1f64);
    let x_ad2 = AD2(2f64, 1f64, 0f64);

    f(x_f64).print();
    f(x_ad1).print();
    f(x_ad2).print();
}

fn f<T: Real>(x: T) -> T {
    return x.powi(2);
}
