extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_print() {
    let mut df = DataFrame::new(vec![]);
    df.push("x", Series::new(vec![1]));
    df.push("y", Series::new(vec![true, false]));
    df.push("z", Series::new(seq(0, 10, 0.01)));

    df.print();
}