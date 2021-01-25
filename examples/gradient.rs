use peroxide::fuga::*;

fn main() {
    f(2f64).print();        // x^3     = 8
    f_grad(2f64).print();   // 3 * x^2 = 12
    f_hess(2f64).print();   // 6 * x   = 12
}

#[ad_function]              // generates f_grad, f_hess
fn f(x: f64) -> f64 {
    x.powi(3)               // x^3
}
