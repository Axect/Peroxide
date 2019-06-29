extern crate peroxide;
use peroxide::*;
use peroxide::operation::union::Number::{F, D};

#[test]
fn test_ops_btw_union_f() {
    let x_f = F(3f64);
    let y_f = F(2f64);

    assert_eq!(x_f + y_f, F(5f64));
    assert_eq!(x_f - y_f, F(1f64));
    assert_eq!(x_f * y_f, F(6f64));
}

#[test]
fn test_ops_btw_union_d() {
    let x_d = D(dual(1, 1)); // y=x
    let y_d = D(dual(1, 2)); // y=x^2

    assert_eq!(x_d + y_d, D(dual(2, 3)));
    assert_eq!(x_d - y_d, D(dual(0, -1)));
    assert_eq!(x_d * y_d, D(dual(1, 3)));
}
