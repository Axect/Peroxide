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

#[cfg(feature = "complex")]
#[test]
fn test_trace_hermitian_real_imag() {
    use num_complex::Complex64;
    use peroxide::complex::matrix::cmatrix;

    let a = cmatrix(
        vec![
            Complex64::new(1f64, 1f64),
            Complex64::new(2f64, -1f64),
            Complex64::new(0f64, 2f64),
            Complex64::new(4f64, -1f64),
        ],
        2,
        2,
        Row,
    );

    // trace = (1+i) + (4-i) = 5
    assert_eq!(a.trace(), Complex64::new(5f64, 0f64));

    // hermitian conjugate: transpose + conjugate
    let ah = a.h();
    assert_eq!(ah[(0, 0)], Complex64::new(1f64, -1f64));
    assert_eq!(ah[(0, 1)], Complex64::new(0f64, -2f64));
    assert_eq!(ah[(1, 0)], Complex64::new(2f64, 1f64));
    assert_eq!(ah[(1, 1)], Complex64::new(4f64, 1f64));

    // (A^H)^H == A
    let ahh = ah.h();
    for i in 0..2 {
        for j in 0..2 {
            assert_eq!(ahh[(i, j)], a[(i, j)]);
        }
    }

    // real / imag parts
    let re = a.real();
    let im = a.imag();
    assert_eq!(re[(0, 1)], 2f64);
    assert_eq!(im[(0, 1)], -1f64);
    assert_eq!(re[(1, 0)], 0f64);
    assert_eq!(im[(1, 0)], 2f64);
}
