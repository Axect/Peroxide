extern crate peroxide;
#[allow(unused_imports)]
use peroxide::fuga::*;

#[cfg(feature = "complex")]
#[test]
fn test_seq() {
    use num_complex::Complex64;
    use peroxide::complex::matrix::cmatrix;

    let v1 = cmatrix(
        vec![
            Complex64::new(1f64, 1f64),
            Complex64::new(2f64, 2f64),
            Complex64::new(3f64, 3f64),
            Complex64::new(4f64, 4f64),
        ],
        2,
        2,
        Row,
    );
    assert_eq!(v1.as_slice()[0], Complex64::new(1f64, 1f64));
}
