extern crate peroxide;
use peroxide::fuga::*;

#[test]
pub fn test_apply() {
    let a = ml_matrix("1 2 3;4 5 6");
    let b = c!(1, 2, 3);
    let c = a.apply(&b);

    assert_eq!(c, c!(14, 32));
}
