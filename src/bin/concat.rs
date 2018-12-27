extern crate peroxide;
use peroxide::*;

fn main() {
    let a = c!(1,2,3);
    let b = c!(4,5,6,7);

    concat(a, b).print();

    let d1 = vec![dual(1, 1)];
    let d2 = vec![dual(2, 2)];
    concat(d1, d2).print();
}
