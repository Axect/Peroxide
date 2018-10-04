extern crate peroxide;

use peroxide::ozone::matrix::*;

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