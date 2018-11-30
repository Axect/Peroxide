extern crate peroxide;
use peroxide::*;

fn main() {
    let a = c!(-1,0,1);
    let b = c!(1,0,1);
    let l = lagrange_polynomial(a, b);
    l.print();
    l.eval_vec(seq!(0, 1, 0.1)).print();

    let u = poly(c!(1,2,3));
    let v = poly(c!(1));
    (u * v).print();
}