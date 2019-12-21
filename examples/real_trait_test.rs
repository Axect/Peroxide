extern crate peroxide;
use peroxide::*;

fn main() {
    let x_f64 = 2f64;
    let x_dual = dual(2, 1);
    let x_hyper = hyper_dual(2, 1, 0);

    f(x_f64).print();
    f(x_dual).print();
    f(x_hyper).print();
}

fn f<T: Real>(x: T) -> T {
    return x.powi(2);
}
