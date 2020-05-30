extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn change_shape_test() {
    let a = matrix(vec![1, 2, 3, 4], 2, 2, Row);
    assert_eq!(a.shape, Row);
    let b = a.change_shape();
    assert_eq!(b.shape, Col);
}
