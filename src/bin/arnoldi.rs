extern crate peroxide;
use peroxide::*;

fn main() {
    let v = normalized_vec(2);
    v.data.print();
    v.data.norm().print();
    v.norm(Frobenius).print();
}

fn normalized_vec(n: usize) -> Matrix {
    let v = rand(n, 1);
    v.fmap(|x| x / v.data.norm())
}
