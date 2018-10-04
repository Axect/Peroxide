extern crate peroxide;

use peroxide::*;

#[test]
fn test_add_matrix() {
    let a = Matrix::new(
        (1..101).collect::<Vec<u32>>(),
        10,
        10,
        Row,
    );
    let b = a.fmap(|x| 2.0*x);
    assert_eq!(a.clone() + a, b);
}

#[test]
fn test_add_f64() {
    let a = Matrix::new(
        vec![1,2,3,4],
        2,
        2,
        Row,
    );
    let b = Matrix::new(
        vec![2,3,4,5],
        2,
        2,
        Row,
    );
    assert_eq!(a + 1.0, b);
}