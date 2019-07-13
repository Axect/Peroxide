extern crate peroxide;
use peroxide::*;

#[test]
fn test_seq() {
    let v1 = c!(2, 4, 6, 8);
    let v2 = seq(2, 8, 2);
    assert_eq!(v1, v2);
}

#[test]
fn test_c() {
    let v1: Vec<f64> = vec![1f64, 2f64, 3f64, 4f64];
    let v2 = c!(1, 2, 3, 4);
    assert_eq!(v1, v2);
}

#[test]
fn test_zeros() {
    let v = zeros!(5);
    assert_eq!(v, vec![0f64; 5]);
}

#[test]
fn test_accumulation() {
    let v1 = c!(1, 2, 3, 4);
    let v2 = seq(5, 8, 1);
    assert_eq!(seq(1, 8, 1), c!(v1; v2));
}

#[test]
fn test_add_matrix() {
    let a = matrix!(1;101;1, 10, 10, Row);
    let b = a.fmap(|x| 2.0 * x);
    assert_eq!(a.clone() + a, b);
}

#[test]
fn test_add_f64() {
    let a = matrix(c!(1, 2, 3, 4), 2, 2, Row);
    let b = matrix(c!(2, 3, 4, 5), 2, 2, Row);
    assert_eq!(a + 1.0, b);
}

#[test]
fn test_col() {
    let a = matrix(seq(1, 4, 1), 2, 2, Row);
    assert_eq!(a.col(0), c!(1, 3));
}

#[test]
fn test_row() {
    let a = matrix!(1;4;1, 2, 2, Row);
    assert_eq!(a.row(0), c!(1, 2));
}

#[test]
fn test_print() {
    let op = Bernoulli(0);
    op.print();
}
