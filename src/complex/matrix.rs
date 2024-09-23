use std::fmt;

use num_complex::Complex;

use crate::fuga::Shape;

/// R-like complex matrix structure
///
/// # Examples
///
/// ```rust
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::ComplexMatrix;
///
/// let v1 = ComplexMatrix {
/// data: vec![
///     Complex64::new(1f64, 1f64),
///     Complex64::new(2f64, 2f64),
///     Complex64::new(3f64, 3f64),
///     Complex64::new(4f64, 4f64),
/// ],
/// row: 2,
/// col: 2,
/// shape: Row,
/// }; // [[1+1i,2+2i],[3+3i,4+4i]]
/// ```
#[derive(Debug, Clone, Default)]
pub struct ComplexMatrix {
    pub data: Vec<Complex<f64>>,
    pub row: usize,
    pub col: usize,
    pub shape: Shape,
}

// =============================================================================
// Various complex matrix constructor
// =============================================================================

/// R-like complex matrix constructor
///
/// # Examples
/// ```rust
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::complex_matrix;
///
/// fn main() {
///     let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
///                       Complex64::new(2f64, 2f64),
///                       Complex64::new(3f64, 3f64),
///                       Complex64::new(4f64, 4f64)],
///                    2, 2, Row
///     );
///     a.col.print(); // Print matrix column
/// }
/// ```
pub fn complex_matrix<T>(v: Vec<T>, r: usize, c: usize, shape: Shape) -> ComplexMatrix
where
    T: Into<Complex<f64>>,
{
    ComplexMatrix {
        data: v
            .into_iter()
            .map(|t| t.into())
            .collect::<Vec<Complex<f64>>>(),
        row: r,
        col: c,
        shape,
    }
}

/// R-like complex matrix constructor (Explicit ver.)
pub fn r_complex_matrix<T>(v: Vec<T>, r: usize, c: usize, shape: Shape) -> ComplexMatrix
where
    T: Into<Complex<f64>>,
{
    complex_matrix(v, r, c, shape)
}

/// Python-like complex matrix constructor
///
/// # Examples
/// ```rust
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
///
/// fn main() {
///     let a = py_complex_matrix(vec![vec![Complex64::new(1f64, 1f64),
///                                         Complex64::new(2f64, 2f64)],
///                                    vec![Complex64::new(3f64, 3f64),
///                                         Complex64::new(4f64, 4f64)]
///     ]);
///     let b = complex_matrix(vec![Complex64::new(1f64, 1f64),
///                                 Complex64::new(2f64, 2f64),
///                                 Complex64::new(3f64, 3f64),
///                                 Complex64::new(4f64, 4f64)],
///                             2, 2, Row
///     );
///     assert_eq!(a, b);
/// }
/// ```
pub fn py_complex_matrix<T>(v: Vec<Vec<T>>) -> ComplexMatrix
where
    T: Into<Complex<f64>> + Copy,
{
    let r = v.len();
    let c = v[0].len();
    let data: Vec<T> = v.into_iter().flatten().collect();
    complex_matrix(data, r, c, Shape::Row)
}

/// Matlab-like matrix constructor
///
/// Note that the entries to the `ml_complex_matrix`
/// needs to be in the `a+bi` format
/// without any spaces between the real and imaginary
/// parts of the Complex number.
///
/// # Examples
/// ```rust
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
///
/// fn main() {
///     let a = ml_complex_matrix("1.0+1.0i 2.0+2.0i;
///                                3.0+3.0i 4.0+4.0i");
///     let b = complex_matrix(vec![Complex64::new(1f64, 1f64),
///                                 Complex64::new(2f64, 2f64),
///                                 Complex64::new(3f64, 3f64),
///                                 Complex64::new(4f64, 4f64)],
///                             2, 2, Row
///     );
///     assert_eq!(a, b);
/// }
/// ```
pub fn ml_complex_matrix(s: &str) -> ComplexMatrix {
    let str_row = s.split(";").collect::<Vec<&str>>();
    let r = str_row.len();
    let str_data = str_row
        .iter()
        .map(|x| x.trim().split(" ").collect::<Vec<&str>>())
        .collect::<Vec<Vec<&str>>>();
    let c = str_data[0].len();
    let data = str_data
        .iter()
        .flat_map(|t| {
            t.iter()
                .map(|x| x.parse::<Complex<f64>>().unwrap())
                .collect::<Vec<Complex<f64>>>()
        })
        .collect::<Vec<Complex<f64>>>();

    complex_matrix(data, r, c, Shape::Row)
}

impl ComplexMatrix {}

// ///  Pretty Print
// impl fmt::Display for ComplexMatrix {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result {
//         write!(f, "{}", self.spread());
//     }
// }

/// PartialEq implements
impl PartialEq for ComplexMatrix {
    fn eq(&self, other: &ComplexMatrix) -> bool {
        if self.shape == other.shape {
            return self.data == other.data;
        } else {
            return false; // todo! self.eq(&other.change_shape())
        }
    }
}
