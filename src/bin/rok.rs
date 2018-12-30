extern crate peroxide;
use peroxide::*;

fn main() {
    let init_vec = c!(0, 1, 1);
    let mut new_vec = one_step_rok4a(init_vec, f, 1e-3);
    new_vec.print();

    for _i in 0 .. 10 {
        new_vec = one_step_rok4a(new_vec, f, 1e-3);
        new_vec.print();
    }
}

fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let t = xs[0];
    let x = xs[1];
    let y = xs[2];

    vec![x * t.cos().pow(2), y * t.sin().pow(2)]
}
