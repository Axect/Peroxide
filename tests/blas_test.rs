extern crate peroxide;
use peroxide::*;

#[test]
#[cfg(feature = "native")]
fn daxpy_test() {
    let a = ml_matrix("1 2; 3 4");
    let b = matrix(vec![1, 3, 2, 4], 2, 2, Col);
    let c = ml_matrix("2 4;6 8");

    assert_eq!(&a + &b, c.clone());
    assert_eq!(a + b, c);
}

#[test]
#[cfg(feature = "native")]
fn dgemv_test() {
    let a = ml_matrix("1 2;3 4");
    let b = c![1, 2];
    let c = ml_matrix("5; 11");
    assert_eq!(&a * &b, c.clone());
    assert_eq!(&(a.change_shape()) * &b, c);
}

#[test]
#[cfg(feature = "native")]
fn dgemm_test() {
    // 2x2
    let a = ml_matrix("1 2;3 4");
    let b = ml_matrix("1 2;3 4");

    let rr = &a * &b;
    let rc = &a * &(b.change_shape());
    let cr = &(a.change_shape()) * &b;
    let cc = &(a.change_shape()) * &(b.change_shape());

    assert_eq!(rc.clone(), cr.clone());
    assert_eq!(rr, rc);
    assert_eq!(cr, cc);

    // 3x3
    let a2 = ml_matrix("1 2 3;4 5 6;7 8 9");
    let b2 = ml_matrix("1 2 3;4 5 6;7 8 9");

    let rr2 = &a2 * &b2;
    let rc2 = &a2 * &(b2.change_shape());
    let cr2 = &(a2.change_shape()) * &b2;
    let cc2 = &(a2.change_shape()) * &(b2.change_shape());

    assert_eq!(rc2.clone(), cr2.clone());
    assert_eq!(rr2, rc2);
    assert_eq!(cr2, cc2);

    // 3x2, 2x3
    let a3 = ml_matrix("1 2;3 4;5 6");
    let b3 = ml_matrix("1 2 3;4 5 6");

    let rr3 = &a3 * &b3;
    let rc3 = &a3 * &(b3.change_shape());
    let cr3 = &(a3.change_shape()) * &b3;
    let cc3 = &(a3.change_shape()) * &(b3.change_shape());

    assert_eq!(rc3.clone(), cr3.clone());
    assert_eq!(rr3, rc3);
    assert_eq!(cr3, cc3);
}
