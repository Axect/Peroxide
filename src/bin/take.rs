extern crate peroxide;
use peroxide::*;

fn main() {
    let a: Matrix = MATLAB::new("1 2 3;4 5 6;7 8 9");
    let b = a.take(2, Row);
    b.print();

    let c = a.take(2, Col);
    c.print();

    let p = a.skip(2, Row);
    p.print();

    let q = a.skip(2, Col);
    q.print();
}