extern crate peroxide;
use peroxide::*;

fn main() {
    let x = 0f64;
    let o1 = ODE { f: test_fun };
    (o1.f)(x).print(); // Now, o1.f becomes f64 -> f64

    let f = o1.f;
    f(1f64).print();

    let o2 = ODE { f: test_fun };

    let g = o2.f;
    g(Dual::new(2, 1)).print(); // Now o2.f becomes Dual -> Dual
}

struct ODE<T: Real> {
    f: fn(T) -> T,
}

fn test_fun<T: Real>(x: T) -> T {
    x + 1f64
}