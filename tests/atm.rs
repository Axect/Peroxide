extern crate peroxide;

use peroxide::*;

#[test]
fn test_add_matrix() {
    let a = matrix!(1;101;1, 10, 10, Row);
    let b = a.fmap(|x| 2.0*x);
    assert_eq!(a.clone() + a, b);
}

#[test]
fn test_add_f64() {
    let a = matrix(
        c!(1,2,3,4),
        2,
        2,
        Row,
    );
    let b = matrix(
        c!(2,3,4,5),
        2,
        2,
        Row,
    );
    assert_eq!(a + 1.0, b);
}

#[test]
fn test_col() {
    let a = matrix(
        seq!(1,4,1),
        2,
        2,
        Row,
    );
    assert_eq!(a.col(0), matrix(c!(1,3),2,1,Col));
}

#[test]
fn test_row() {
    let a = matrix!(1;4;1, 2, 2, Row);
    assert_eq!(a.row(0), matrix(c!(1,2), 1, 2, Row));
}