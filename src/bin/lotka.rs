extern crate peroxide;
use peroxide::*;

#[allow(unused_must_use)]
fn main() {
    // t = 0, x = 2, y = 1
    let mut t = 0f64;
    let xs = c!(2, 1);
    let rk_records = solve(lotka_volterra, xs.clone(), (0, 10), 1e-3, RK4);
    let bdf_records = solve(lotka_volterra, xs, (0, 10), 1e-3, BDF1(1e-15));
    rk_records.print();
    bdf_records.print();
}

fn lotka_volterra(t: Dual, xs: Vec<Dual>) -> Vec<Dual> {
    let a = 4.;
    let c = 1.;

    let x = xs[0];
    let y = xs[1];

    vec![
        a * (x - x*y),
        -c * (y - x*y)
    ]
}
