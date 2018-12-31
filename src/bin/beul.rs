extern crate peroxide;
use peroxide::*;

fn main() {
    // t = 0, x = 2, y = 1
    let mut xs = c!(0, 2, 1);
    let mut records = zeros(10001, 3);
    records.subs_row(0, xs.clone());
    
    for i in 1 .. 10001 {
        xs = one_step_bdf1(xs.clone(), lotka_volterra, 1e-3, 1e-15);
        records.subs_row(i, xs.clone());
    }

    records.write("beul.csv");
}

fn lotka_volterra(xs: Vec<Dual>) -> Vec<Dual> {
    let a = 4.;
    let c = 1.;

    let x = xs[1];
    let y = xs[2];

    vec![
        a * (x.clone() - x.clone() * y.clone()),
        -c * (y.clone() - x * y)
    ]
}
