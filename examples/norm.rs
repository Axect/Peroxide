extern crate peroxide;
use peroxide::*;
use peroxide::traits::math::{Norm, Normed};

fn main() {
    let a = matrix(c!(1, 1, 1, 1), 2, 2, Row);
    a.norm(Norm::F).print();
    a.norm(Norm::L1).print();
    a.norm(Norm::LInf).print();
    a.norm(Norm::Lpq(2f64, 2f64)).print();
    a.norm(Norm::Lpq(2f64, 1f64)).print();
}
