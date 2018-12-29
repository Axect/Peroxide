extern crate peroxide;
use peroxide::*;

fn main() {
    let init_vec = c!(0, 1, 2);
    let (h1, v1, w1) = modified_arnoldi(init_vec.clone(), f);
    h1.print();
    v1.print();
    w1.print();
    let mut new_vec = one_step_rok4a(init_vec, f, 1e-1);

    let (h2, v2, w2) = modified_arnoldi(new_vec.clone(), f);
    h2.print();
    v2.print();
    w2.print();


    let mut second_vec = one_step_rok4a(new_vec, f, 1e-1);
}

fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let t = xs[0];
    let x = xs[1];
    let y = xs[2];

    vec![x* t.cos().pow(2), y * t.sin().pow(2)]
}
