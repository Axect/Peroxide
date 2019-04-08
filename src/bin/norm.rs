extern crate peroxide;
use peroxide::*;

fn main() {
    let a = matrix(c!(1,1,1,1), 2, 2, Row);
    a.norm(Frobenius).print();
    a.norm(One).print();
    a.norm(Infinity).print();
    a.norm(PQ(2,2)).print();
    a.norm(PQ(2,1)).print();
}
