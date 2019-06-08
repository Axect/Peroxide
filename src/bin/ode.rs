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

    o1.eval_f64(1f64).print();
    o1.eval_dual(Dual::new(2, 1)).print();
    o2.eval_f64(1f64).print();
    o2.eval_dual(Dual::new(2, 1)).print();
}

struct ODE<T: Real> {
    f: fn(T) -> T,
}

fn test_fun<T: Real>(x: T) -> T {
    x + 1f64
}

impl<T: Real> ODE<T> {
    fn eval_f64(&self, x: f64) -> f64 {
        (self.f)(T::from_f64(x)).to_f64()
    }

    fn eval_dual(&self, x: Dual) -> Dual {
        (self.f)(T::from_dual(x)).to_dual()
    }
}
