extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_binomial() {
    let b = Binomial(100, 0.8);
    b.sample(10).print();
    assert!(nearly_eq(b.mean(), 80f64));
    assert!(nearly_eq(b.var(), 16f64));
}
