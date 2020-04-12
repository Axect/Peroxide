extern crate peroxide;
use peroxide::*;

#[test]
fn test_qr() {
    let a = ml_matrix("12 -51 4;6 167 -68;-4 24 -41");
    let qr = a.qr();
    assert_eq!(a, qr.q() * qr.r());
}