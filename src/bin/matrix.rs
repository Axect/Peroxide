extern crate peroxide;
use peroxide::*;

fn main() {
    let p = py_matrix(vec![c!(1, 2), c!(3, 4)]);
    p.print();

    let m = ml_matrix("1 2 3; 3 4 5; 6 7 8");
    m.print();
}