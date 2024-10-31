extern crate peroxide;
#[allow(unused_imports)]
use peroxide::fuga::*;

#[cfg(feature = "complex")]
#[test]
fn test_seq() {
    use num_complex::Complex64;
    use peroxide::complex::matrix::ComplexMatrix;

    let v1 = ComplexMatrix {
        data: vec![
            Complex64::new(1f64, 1f64),
            Complex64::new(2f64, 2f64),
            Complex64::new(3f64, 3f64),
            Complex64::new(4f64, 4f64),
        ],
        row: 2,
        col: 2,
        shape: Row,
    };
    assert_eq!(v1.data[0], Complex64::new(1f64, 1f64));
}
