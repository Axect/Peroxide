extern crate peroxide;
use peroxide::*;

fn main() {
    let a: Matrix = MATLAB::new("1 2; 3 4");
    a.print();

    let b: Matrix = PYTHON::new(vec![c!(1, 2), c!(3, 4)]);
    b.print();

    let c: Matrix = R::new(c!(1,2,3,4), 2, 2, Row);
    c.print();
}