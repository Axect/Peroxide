extern crate peroxide;
use peroxide::*;

fn main() {
    let a = matrix(c!(1,2,3,4,5,6), 3, 2, Row);
    (&a * &a.t()).print();
}