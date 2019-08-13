#[macro_use]
extern crate peroxide;
use peroxide::operation::number::Number::{D, F};
use peroxide::*;

#[test]
fn test_ops_btw_number_f() {
    let x_f = F(3f64);
    let y_f = F(2f64);

    assert_eq!(x_f + y_f, F(5f64));
    assert_eq!(x_f - y_f, F(1f64));
    assert_eq!(x_f * y_f, F(6f64));
}

#[test]
fn test_ops_btw_number_d() {
    let x_d = D(dual(1, 1)); // y=x
    let y_d = D(dual(1, 2)); // y=x^2

    assert_eq!(x_d + y_d, D(dual(2, 3)));
    assert_eq!(x_d - y_d, D(dual(0, -1)));
    assert_eq!(x_d * y_d, D(dual(1, 3)));
}

#[test]
fn test_function_for_number() {
    let fl = F(2f64);
    let _du = D(dual(1, 1));

    assert_eq!(f(fl), F(14f64));
}

fn f(n: Number) -> Number {
    2f64 * n * 3f64 + 2f64
}
