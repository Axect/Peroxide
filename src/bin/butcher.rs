extern crate peroxide;
use peroxide::*;
use std::ops::Add;

fn main() {
    // RK4 Matrix
    let mut a = zeros(4, 4);
    a[(1,0)] = 0.5;
    a[(2,1)] = 0.5;
    a[(3,2)] = 1f64;
    a.print();

    let c = c!(0, 0.5, 0.5, 1);
    let b = c!(1f64/6f64, 1f64/3f64, 1f64/3f64, 1f64/6f64);

    b.print();
    c.print();
}

fn f<T: Real>(t: T, x: T) -> T {
    t + x
}