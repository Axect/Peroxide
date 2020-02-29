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
    let a = matrix!(1;100;1, 10, 10, Row);
    assert_eq!(&a + &a, 2f64 * a);
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

#[test]
fn test_row_map() {
    let m = ml_matrix("1 2;3 4");
    let n = m.row_map(|v| v.normalize());
    let o = matrix(
        vec![
            1f64 / 5f64.sqrt(),
            2f64 / 5f64.sqrt(),
            3f64 / 5f64,
            4f64 / 5f64,
        ],
        2,
        2,
        Row,
    );
    assert_eq!(n, o);
}

#[test]
fn test_col_map() {
    let m = ml_matrix("1 3;2 4");
    let n = m.col_map(|v| v.normalize());
    let o = matrix(
        vec![
            1f64 / 5f64.sqrt(),
            2f64 / 5f64.sqrt(),
            3f64 / 5f64,
            4f64 / 5f64,
        ],
        2,
        2,
        Col,
    );
    assert_eq!(n, o);
}
