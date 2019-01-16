extern crate peroxide;
use peroxide::*;

#[allow(unused_must_use)]
fn main() {
    // t = 0, x = 2, y = 1
    let mut t = 0f64;
    let mut xs = c!(2, 1);
    let rk_records = rk4(t, xs.clone(), lotka_volterra, 1e-3, 10000);
    rk_records.write_with_header("example_data/lotka_rk.csv", vec!["t", "x", "y"]);
}

fn lotka_volterra(t: f64, xs: Vec<f64>) -> Vec<f64> {
    let a = 4.;
    let c = 1.;

    let x = xs[0];
    let y = xs[1];

    vec![
        a * (x - x * y),
        -c * (y - x * y)
    ]
}


fn lotka_volterra_dual(xs: Vec<Dual>) -> Vec<Dual> {
    let a = 4.;
    let c = 1.;

    let x = xs[1];
    let y = xs[2];

    vec![
        a * (x.clone() - x.clone() * y.clone()),
        -c * (y.clone() - x * y)
    ]
}
