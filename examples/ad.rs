extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = AD1::new(2f64, 1f64);
    a.print();
    let a2 = AD2::from(a);

    let b = AD2::new(4f64, 4f64, 2f64);
    b.print();

    (a + b).print();
    (a - b).print();
    (a * b).print();
    (a / b).print();
    (a2 / b).print();

    let b_iter = b.iter();
    println!("{:?}", b_iter.size_hint());
    let b_skip = b_iter.skip(1);
    println!("{:?}", b_skip.size_hint());
    let b_take = b_skip.take(2);
    println!("{:?}", b_take.size_hint());
    println!("{:?}", b_take.rev().size_hint());

    let b_iter = b.iter();
    println!("{:?}", b_iter.len());
    let b_skip = b_iter.skip(1);
    println!("{:?}", b_skip.len());
    let b_take = b_skip.take(2);
    println!("{:?}", b_take.len());
    println!("{:?}", b_take.rev().len());

    a.ln().print();
    b.ln().print();
}
