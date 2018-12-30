extern crate peroxide;
use peroxide::*;

fn main() {
    // t = 0, x = 2, y = 1
    let mut xs = c!(0, 2, 1);
    
    for _i in 0 .. 100 {
        xs = one_step_bdf1(xs.clone(), lotka_volterra, 1e-1, 1e-15);
        xs.print();
    }
}

fn lotka_volterra(xs: Vec<Dual>) -> Vec<Dual> {
    let a = 4.;
    let c = 1.;

    let t = xs[0];
    let x = xs[1];
    let y = xs[2];

    vec![
        a * (x.clone() - x.clone() * y.clone()),
        -c * (y.clone() - x * y)
    ]
}