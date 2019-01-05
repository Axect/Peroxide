extern crate peroxide;
use peroxide::*;

#[allow(unused_must_use)]
fn main() {
    // t = 0, x = 2, y = 1
    let mut xs = c!(0, 2, 1);
    let bdf_records = bdf1(xs.clone(), lotka_volterra_dual, 1e-3, 1e-15, 10000);
    let rk_records = rk4(xs.clone(), lotka_volterra, 1e-3, 10000);

    bdf_records.write("example_data/lotka_bdf.csv");
    rk_records.write("example_data/lotka_rk.csv");
}

fn lotka_volterra(xs: Vec<f64>) -> Vec<f64> {
    let a = 4.;
    let c = 1.;

    let x = xs[1];
    let y = xs[2];

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
