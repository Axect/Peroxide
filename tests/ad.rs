extern crate peroxide;
#[allow(unused_imports)]
use peroxide::fuga::*;

#[test]
fn add_test() {
    let a = AD::from([3f64, 1f64, 0f64]);
    let b = AD::from([2f64, 1f64, 0f64]);
    assert_eq!(&a + &b, AD::from([5f64, 2f64, 0f64]));
}