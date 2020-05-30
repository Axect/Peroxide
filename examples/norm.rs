extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = ml_matrix("1 1; 1 1");
    a.norm(Norm::F).print();
    a.norm(Norm::L1).print();
    a.norm(Norm::LInf).print();
    a.norm(Norm::Lpq(2f64, 2f64)).print();
    a.norm(Norm::Lpq(2f64, 1f64)).print();
}
