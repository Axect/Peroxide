extern crate peroxide;
use peroxide::*;

fn main() {
    let init_vec = c!(1,1);
    let (h, v, w) = modified_arnoldi(init_vec, f);
    h.print();
    v.print();
    w.print();
}

fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let t = xs[0];
    let x = xs[1];

    vec![t.pow(2) * 3. * x]
}
