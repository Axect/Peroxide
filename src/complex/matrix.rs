use std::{
    cmp::{max, min},
    fmt,
    ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub},
};

use anyhow::{bail, Result};
use matrixmultiply::CGemmOption;
use num_complex::Complex;
use rand_distr::num_traits::{One, Zero};

use crate::{
    fuga::{
        nearly_eq, tab, ConcatenateError, InnerProduct, LinearOp, MatrixProduct, Norm, Normed,
        Shape, Vector,
    },
    traits::{fp::FPMatrix, mutable::MutMatrix},
};

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

///  Pretty Print
impl fmt::Display for ComplexMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

/// PartialEq implements
impl PartialEq for ComplexMatrix {
    fn eq(&self, other: &ComplexMatrix) -> bool {
        if self.shape == other.shape {
            self.data
                .clone()
                .into_iter()
                .zip(other.data.clone())
                .all(|(x, y)| nearly_eq(x.re, y.re) && nearly_eq(x.im, y.im))
                && self.row == other.row
        } else {
            self.eq(&other.change_shape())
        }
    }
}

impl ComplexMatrix {
    /// Raw pointer for `self.data`
    pub fn ptr(&self) -> *const Complex<f64> {
        &self.data[0] as *const Complex<f64>
    }

    /// Raw mutable pointer for `self.data`
    pub fn mut_ptr(&mut self) -> *mut Complex<f64> {
        &mut self.data[0] as *mut Complex<f64>
    }

    /// Slice of `self.data`
    ///
    /// # Examples
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    ///
    /// let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                             Complex64::new(2f64, 2f64),
    ///                             Complex64::new(3f64, 3f64),
    ///                             Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    /// let b = a.as_slice();
    /// assert_eq!(b, &[Complex64::new(1f64, 1f64),
    ///                 Complex64::new(2f64, 2f64),
    ///                 Complex64::new(3f64, 3f64),
    ///                 Complex64::new(4f64, 4f64)]);
    /// ```
    pub fn as_slice(&self) -> &[Complex<f64>] {
        &self.data[..]
    }

    /// Mutable slice of `self.data`
    ///
    /// # Examples
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    ///
    /// let mut a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    /// let mut b = a.as_mut_slice();
    /// b[1] = Complex64::new(5f64, 5f64);
    /// assert_eq!(b, &[Complex64::new(1f64, 1f64),
    ///                 Complex64::new(5f64, 5f64),
    ///                 Complex64::new(3f64, 3f64),
    ///                 Complex64::new(4f64, 4f64)]);
    /// assert_eq!(a, complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                   Complex64::new(5f64, 5f64),
    ///                                   Complex64::new(3f64, 3f64),
    ///                                   Complex64::new(4f64, 4f64)],
    ///                               2, 2, Row));
    /// ```
    pub fn as_mut_slice(&mut self) -> &mut [Complex<f64>] {
        &mut self.data[..]
    }

    /// Change Bindings
    ///
    /// `Row` -> `Col` or `Col` -> `Row`
    ///
    /// # Examples
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    ///
    /// let mut a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    /// assert_eq!(a.shape, Row);
    /// let b = a.change_shape();
    /// assert_eq!(b.shape, Col);
    /// ```
    pub fn change_shape(&self) -> Self {
        let r = self.row;
        let c = self.col;
        assert_eq!(r * c, self.data.len());
        let l = r * c - 1;
        let mut data: Vec<Complex<f64>> = self.data.clone();
        let ref_data = &self.data;

        match self.shape {
            Shape::Row => {
                for i in 0..l {
                    let s = (i * c) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                complex_matrix(data, r, c, Shape::Col)
            }
            Shape::Col => {
                for i in 0..l {
                    let s = (i * r) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                complex_matrix(data, r, c, Shape::Row)
            }
        }
    }

    /// Change Bindings Mutably
    ///
    /// `Row` -> `Col` or `Col` -> `Row`
    ///
    /// # Examples
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    ///
    /// let mut a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    /// assert_eq!(a.shape, Row);
    /// let b = a.change_shape_mut();
    /// assert_eq!(b.shape, Col);
    /// ```
    pub fn change_shape_mut(&mut self) -> Self {
        let r = self.row;
        let c = self.col;
        assert_eq!(r * c, self.data.len());
        let l = r * c - 1;
        let mut data: Vec<Complex<f64>> = self.data.clone();
        let ref_data = &self.data;

        match self.shape {
            Shape::Row => {
                for i in 0..l {
                    let s = (i * c) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                complex_matrix(data, r, c, Shape::Col)
            }
            Shape::Col => {
                for i in 0..l {
                    let s = (i * r) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                complex_matrix(data, r, c, Shape::Row)
            }
        }
    }

    /// Spread data(1D vector) to 2D formatted String
    ///
    /// # Examples
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    ///
    /// let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    /// println!("{}", a.spread()); // same as println!("{}", a);
    /// // Result:
    /// //       c[0]    c[1]
    /// // r[0]  1+1i    3+3i
    /// // r[1]  2+2i    4+4i
    /// ```
    pub fn spread(&self) -> String {
        assert_eq!(self.row * self.col, self.data.len());
        let r = self.row;
        let c = self.col;
        let mut key_row = 20usize;
        let mut key_col = 20usize;

        if r > 100 || c > 100 || (r > 20 && c > 20) {
            let part = if r <= 10 {
                key_row = r;
                key_col = 100;
                self.take_col(100)
            } else if c <= 10 {
                key_row = 100;
                key_col = c;
                self.take_row(100)
            } else {
                self.take_row(20).take_col(20)
            };
            return format!(
                "Result is too large to print - {}x{}\n only print {}x{} parts:\n{}",
                self.row.to_string(),
                self.col.to_string(),
                key_row.to_string(),
                key_col.to_string(),
                part.spread()
            );
        }

        // Find maximum length of data
        let sample = self.data.clone();
        let mut space: usize = sample
            .into_iter()
            .map(
                |x| min(format!("{:.4}", x).len(), x.to_string().len()), // Choose minimum of approx vs normal
            )
            .fold(0, |x, y| max(x, y))
            + 1;

        if space < 5 {
            space = 5;
        }

        let mut result = String::new();

        result.push_str(&tab("", 5));
        for i in 0..c {
            result.push_str(&tab(&format!("c[{}]", i), space)); // Header
        }
        result.push('\n');

        for i in 0..r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for j in 0..c {
                let st1 = format!("{:.4}", self[(i, j)]); // Round at fourth position
                let st2 = self[(i, j)].to_string(); // Normal string
                let mut st = st2.clone();

                // Select more small thing
                if st1.len() < st2.len() {
                    st = st1;
                }

                result.push_str(&tab(&st, space));
            }
            if i == (r - 1) {
                break;
            }
            result.push('\n');
        }

        return result;
    }

    /// Extract Column
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
    ///     let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                             Complex64::new(2f64, 2f64),
    ///                             Complex64::new(3f64, 3f64),
    ///                             Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///         );
    ///     assert_eq!(a.col(0), vec![Complex64::new(1f64, 1f64), Complex64::new(3f64, 3f64)]);
    /// }
    /// ```
    pub fn col(&self, index: usize) -> Vec<Complex<f64>> {
        assert!(index < self.col);
        let mut container: Vec<Complex<f64>> = vec![Complex::zero(); self.row];
        for i in 0..self.row {
            container[i] = self[(i, index)];
        }
        container
    }

    /// Extract Row
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
    ///     let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                             Complex64::new(2f64, 2f64),
    ///                             Complex64::new(3f64, 3f64),
    ///                             Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///         );
    ///     assert_eq!(a.row(0), vec![Complex64::new(1f64, 1f64), Complex64::new(2f64, 2f64)]);
    /// }
    /// ```
    pub fn row(&self, index: usize) -> Vec<Complex<f64>> {
        assert!(index < self.row);
        let mut container: Vec<Complex<f64>> = vec![Complex::zero(); self.col];
        for i in 0..self.col {
            container[i] = self[(index, i)];
        }
        container
    }

    /// Extract diagonal components
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
    ///     let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///          );
    ///     assert_eq!(a.diag(), vec![Complex64::new(1f64, 1f64) ,Complex64::new(4f64, 4f64)]);
    /// }
    /// ```
    pub fn diag(&self) -> Vec<Complex<f64>> {
        let mut container = vec![Complex::zero(); self.row];
        let r = self.row;
        let c = self.col;
        assert_eq!(r, c);

        let c2 = c + 1;
        for i in 0..r {
            container[i] = self.data[i * c2];
        }
        container
    }

    /// Transpose
    ///
    /// # Examples
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    ///
    /// let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                             Complex64::new(2f64, 2f64),
    ///                             Complex64::new(3f64, 3f64),
    ///                             Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    /// let a_t = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                               Complex64::new(2f64, 2f64),
    ///                               Complex64::new(3f64, 3f64),
    ///                               Complex64::new(4f64, 4f64)],
    ///                             2, 2, Col
    ///     );
    ///
    /// assert_eq!(a.transpose(), a_t);
    /// ```
    pub fn transpose(&self) -> Self {
        match self.shape {
            Shape::Row => complex_matrix(self.data.clone(), self.col, self.row, Shape::Col),
            Shape::Col => complex_matrix(self.data.clone(), self.col, self.row, Shape::Row),
        }
    }

    /// R-like transpose function
    ///
    /// # Examples
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    ///
    /// let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                             Complex64::new(2f64, 2f64),
    ///                             Complex64::new(3f64, 3f64),
    ///                             Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    /// let a_t = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                               Complex64::new(2f64, 2f64),
    ///                               Complex64::new(3f64, 3f64),
    ///                               Complex64::new(4f64, 4f64)],
    ///                             2, 2, Col
    ///     );
    ///
    /// assert_eq!(a.t(), a_t);
    /// ```
    pub fn t(&self) -> Self {
        self.transpose()
    }

    /// Should check shape
    pub fn subs(&mut self, idx: usize, v: &Vec<Complex<f64>>) {
        let p = &mut self.mut_ptr();
        match self.shape {
            Shape::Row => {
                let c = self.col;
                unsafe {
                    p.add(idx * c).copy_from(v.as_ptr(), c);
                }
            }
            Shape::Col => {
                let r = self.row;
                unsafe {
                    p.add(idx * r).copy_from(v.as_ptr(), r);
                }
            }
        }
    }

    /// Substitute Col
    #[inline]
    pub fn subs_col(&mut self, idx: usize, v: &Vec<Complex<f64>>) {
        for i in 0..self.row {
            self[(i, idx)] = v[i];
        }
    }

    /// Substitute Row
    #[inline]
    pub fn subs_row(&mut self, idx: usize, v: &Vec<Complex<f64>>) {
        for j in 0..self.col {
            self[(idx, j)] = v[j];
        }
    }

    /// From index operations
    pub fn from_index<F, G>(f: F, size: (usize, usize)) -> ComplexMatrix
    where
        F: Fn(usize, usize) -> G + Copy,
        G: Into<Complex<f64>>,
    {
        let row = size.0;
        let col = size.1;

        let mut mat = complex_matrix(vec![Complex::zero(); row * col], row, col, Shape::Row);

        for i in 0..row {
            for j in 0..col {
                mat[(i, j)] = f(i, j).into();
            }
        }
        mat
    }

    /// Matrix to `Vec<Vec<Complex<f64>>>`
    ///
    /// To send `Matrix` to `inline-python`
    pub fn to_vec(&self) -> Vec<Vec<Complex<f64>>> {
        let mut result = vec![vec![Complex::zero(); self.col]; self.row];
        for i in 0..self.row {
            result[i] = self.row(i);
        }
        result
    }

    pub fn to_diag(&self) -> ComplexMatrix {
        assert_eq!(self.row, self.col, "Should be square matrix");
        let mut result = complex_matrix(
            vec![Complex::zero(); self.row * self.col],
            self.row,
            self.col,
            Shape::Row,
        );
        let diag = self.diag();
        for i in 0..self.row {
            result[(i, i)] = diag[i];
        }
        result
    }

    /// Submatrix
    ///
    /// # Description
    /// Return below elements of complex matrix to a new complex matrix
    ///
    /// $$
    /// \begin{pmatrix}
    /// \\ddots & & & & \\\\
    ///   & start & \\cdots & end.1 & \\\\
    ///   & \\vdots & \\ddots & \\vdots & \\\\
    ///   & end.0 & \\cdots & end & \\\\
    ///   & & & & \\ddots
    /// \end{pmatrix}
    /// $$
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
    ///     let a = ml_complex_matrix("1.0+1.0i 2.0+2.0i 3.0+3.0i;
    ///                                4.0+4.0i 5.0+5.0i 6.0+6.0i;
    ///                                7.0+7.0i 8.0+8.0i 9.0+9.0i");
    ///     let b = complex_matrix(vec![Complex64::new(5f64, 5f64),
    ///                                 Complex64::new(6f64, 6f64),
    ///                                 Complex64::new(8f64, 8f64),
    ///                                 Complex64::new(9f64, 9f64)],
    ///                             2, 2, Row
    ///     );
    ///     let c = a.submat((1, 1), (2, 2));
    ///     assert_eq!(b, c);
    /// }
    /// ```
    pub fn submat(&self, start: (usize, usize), end: (usize, usize)) -> ComplexMatrix {
        let row = end.0 - start.0 + 1;
        let col = end.1 - start.1 + 1;
        let mut result = complex_matrix(vec![Complex::zero(); row * col], row, col, self.shape);
        for i in 0..row {
            for j in 0..col {
                result[(i, j)] = self[(start.0 + i, start.1 + j)];
            }
        }
        result
    }

    /// Substitute complex matrix to specific position
    ///
    /// # Description
    /// Substitute below elements of complex matrix
    ///
    /// $$
    /// \begin{pmatrix}
    /// \\ddots & & & & \\\\
    ///   & start & \\cdots & end.1 & \\\\
    ///   & \\vdots & \\ddots & \\vdots & \\\\
    ///   & end.0 & \\cdots & end & \\\\
    ///   & & & & \\ddots
    /// \end{pmatrix}
    /// $$
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    ///
    /// fn main() {
    ///     let mut a = ml_complex_matrix("1.0+1.0i 2.0+2.0i 3.0+3.0i;
    ///                                4.0+4.0i 5.0+5.0i 6.0+6.0i;
    ///                                7.0+7.0i 8.0+8.0i 9.0+9.0i");
    ///     let b = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row);
    ///     let c = ml_complex_matrix("1.0+1.0i 2.0+2.0i 3.0+3.0i;
    ///                                4.0+4.0i 1.0+1.0i 2.0+2.0i;
    ///                                7.0+7.0i 3.0+3.0i 4.0+4.0i");
    ///     a.subs_mat((1,1), (2,2), &b);
    ///     assert_eq!(a, c);       
    /// }
    /// ```
    pub fn subs_mat(&mut self, start: (usize, usize), end: (usize, usize), m: &ComplexMatrix) {
        let row = end.0 - start.0 + 1;
        let col = end.1 - start.1 + 1;
        for i in 0..row {
            for j in 0..col {
                self[(start.0 + i, start.1 + j)] = m[(i, j)];
            }
        }
    }
}

// =============================================================================
// Mathematics for Matrix
// =============================================================================
impl Vector for ComplexMatrix {
    type Scalar = Complex<f64>;

    fn add_vec(&self, other: &Self) -> Self {
        assert_eq!(self.row, other.row);
        assert_eq!(self.col, other.col);

        let mut result = complex_matrix(self.data.clone(), self.row, self.col, self.shape);
        for i in 0..self.row {
            for j in 0..self.col {
                result[(i, j)] += other[(i, j)];
            }
        }
        result
    }

    fn sub_vec(&self, other: &Self) -> Self {
        assert_eq!(self.row, other.row);
        assert_eq!(self.col, other.col);

        let mut result = complex_matrix(self.data.clone(), self.row, self.col, self.shape);
        for i in 0..self.row {
            for j in 0..self.col {
                result[(i, j)] -= other[(i, j)];
            }
        }
        result
    }

    fn mul_scalar(&self, other: Self::Scalar) -> Self {
        let scalar = other;
        self.fmap(|x| x * scalar)
    }
}

impl Normed for ComplexMatrix {
    type UnsignedScalar = f64;

    fn norm(&self, kind: Norm) -> Self::UnsignedScalar {
        match kind {
            Norm::F => {
                let mut s = Complex::zero();
                for i in 0..self.data.len() {
                    s += self.data[i].powi(2);
                }
                s.sqrt().re
            }
            Norm::Lpq(p, q) => {
                let mut s = Complex::zero();
                for j in 0..self.col {
                    let mut s_row = Complex::zero();
                    for i in 0..self.row {
                        s_row += self[(i, j)].powi(p as i32);
                    }
                    s += s_row.powf(q / p);
                }
                s.powf(1f64 / q).re
            }
            Norm::L1 => {
                let mut m = Complex::zero();
                match self.shape {
                    Shape::Row => self.change_shape().norm(Norm::L1),
                    Shape::Col => {
                        for c in 0..self.col {
                            let s: Complex<f64> = self.col(c).iter().sum();
                            if s.re > m.re {
                                m = s;
                            }
                        }
                        m.re
                    }
                }
            }
            Norm::LInf => {
                let mut m = Complex::zero();
                match self.shape {
                    Shape::Col => self.change_shape().norm(Norm::LInf),
                    Shape::Row => {
                        for r in 0..self.row {
                            let s: Complex<f64> = self.row(r).iter().sum();
                            if s.re > m.re {
                                m = s;
                            }
                        }
                        m.re
                    }
                }
            }
            Norm::L2 => {
                unimplemented!()
            }
            Norm::Lp(_) => unimplemented!(),
        }
    }

    fn normalize(&self, _kind: Norm) -> Self
    where
        Self: Sized,
    {
        unimplemented!()
    }
}

/// Frobenius inner product
impl InnerProduct for ComplexMatrix {
    fn dot(&self, rhs: &Self) -> Complex<f64> {
        if self.shape == rhs.shape {
            self.data.dot(&rhs.data)
        } else {
            self.data.dot(&rhs.change_shape().data)
        }
    }
}

/// TODO: Transpose

/// Matrix as Linear operator for Vector
#[allow(non_snake_case)]
impl LinearOp<Vec<Complex<f64>>, Vec<Complex<f64>>> for ComplexMatrix {
    fn apply(&self, other: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        assert_eq!(self.col, other.len());
        let mut c = vec![Complex::zero(); self.row];
        complex_gemv(Complex::one(), self, other, Complex::zero(), &mut c);
        c
    }
}

/// R like cbind - concatenate two comlex matrix by column direction
pub fn complex_cbind(m1: ComplexMatrix, m2: ComplexMatrix) -> Result<ComplexMatrix> {
    let mut temp = m1;
    if temp.shape != Shape::Col {
        temp = temp.change_shape();
    }

    let mut temp2 = m2;
    if temp2.shape != Shape::Col {
        temp2 = temp2.change_shape();
    }

    let mut v = temp.data;
    let mut c = temp.col;
    let r = temp.row;

    if r != temp2.row {
        bail!(ConcatenateError::DifferentLength);
    }
    v.extend_from_slice(&temp2.data[..]);
    c += temp2.col;
    Ok(complex_matrix(v, r, c, Shape::Col))
}

/// R like rbind - concatenate two complex matrix by row direction
/// ```
pub fn complex_rbind(m1: ComplexMatrix, m2: ComplexMatrix) -> Result<ComplexMatrix> {
    let mut temp = m1;
    if temp.shape != Shape::Row {
        temp = temp.change_shape();
    }

    let mut temp2 = m2;
    if temp2.shape != Shape::Row {
        temp2 = temp2.change_shape();
    }

    let mut v = temp.data;
    let c = temp.col;
    let mut r = temp.row;

    if c != temp2.col {
        bail!(ConcatenateError::DifferentLength);
    }
    v.extend_from_slice(&temp2.data[..]);
    r += temp2.row;
    Ok(complex_matrix(v, r, c, Shape::Row))
}

impl MatrixProduct for ComplexMatrix {
    fn kronecker(&self, other: &Self) -> Self {
        let r1 = self.row;
        let c1 = self.col;

        let mut result = self[(0, 0)] * other;

        for j in 1..c1 {
            let n = self[(0, j)] * other;
            result = complex_cbind(result, n).unwrap();
        }

        for i in 1..r1 {
            let mut m = self[(i, 0)] * other;
            for j in 1..c1 {
                let n = self[(i, j)] * other;
                m = complex_cbind(m, n).unwrap();
            }
            result = complex_rbind(result, m).unwrap();
        }
        result
    }

    fn hadamard(&self, other: &Self) -> Self {
        assert_eq!(self.row, other.row);
        assert_eq!(self.col, other.col);

        let r = self.row;
        let c = self.col;

        let mut m = complex_matrix(vec![Complex::zero(); r * c], r, c, self.shape);
        for i in 0..r {
            for j in 0..c {
                m[(i, j)] = self[(i, j)] * other[(i, j)]
            }
        }
        m
    }
}

// =============================================================================
// Common Properties of Matrix & Vec<f64>
// =============================================================================
/// `Complex Matrix` to `Vec<Complex<f64>>`
impl Into<Vec<Complex<f64>>> for ComplexMatrix {
    fn into(self) -> Vec<Complex<f64>> {
        self.data
    }
}

/// `&ComplexMatrix` to `&Vec<Complex<f64>>`
impl<'a> Into<&'a Vec<Complex<f64>>> for &'a ComplexMatrix {
    fn into(self) -> &'a Vec<Complex<f64>> {
        &self.data
    }
}

/// `Vec<Complex<f64>>` to `ComplexMatrix`
impl Into<ComplexMatrix> for Vec<Complex<f64>> {
    fn into(self) -> ComplexMatrix {
        let l = self.len();
        complex_matrix(self, l, 1, Shape::Col)
    }
}

/// `&Vec<Complex<f64>>` to `ComplexMatrix`
impl Into<ComplexMatrix> for &Vec<Complex<f64>> {
    fn into(self) -> ComplexMatrix {
        let l = self.len();
        complex_matrix(self.clone(), l, 1, Shape::Col)
    }
}

// =============================================================================
// Standard Operation for Complex Matrix (ADD)
// =============================================================================

/// Element-wise addition of Complex Matrix
///
/// # Caution
/// > You should remember ownership.
/// > If you use ComplexMatrix `a,b` then you can't use them after.
impl Add<ComplexMatrix> for ComplexMatrix {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        assert_eq!(&self.row, &other.row);
        assert_eq!(&self.col, &other.col);

        let mut result = complex_matrix(self.data.clone(), self.row, self.col, self.shape);
        for i in 0..self.row {
            for j in 0..self.col {
                result[(i, j)] += other[(i, j)];
            }
        }
        result
    }
}

impl<'a, 'b> Add<&'b ComplexMatrix> for &'a ComplexMatrix {
    type Output = ComplexMatrix;

    fn add(self, rhs: &'b ComplexMatrix) -> Self::Output {
        self.add_vec(rhs)
    }
}

/// Element-wise addition between Complex Matrix & Complex<f64>
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
///     let mut a = ml_complex_matrix("1.0+1.0i 2.0+2.0i;
///                                    4.0+4.0i 5.0+5.0i");
///     let a_exp = ml_complex_matrix("2.0+2.0i 3.0+3.0i;
///                                    5.0+5.0i 6.0+6.0i");
///     assert_eq!(a + Complex64::new(1_f64, 1_f64), a_exp);
/// }
/// ```
impl<T> Add<T> for ComplexMatrix
where
    T: Into<Complex<f64>> + Copy,
{
    type Output = Self;
    fn add(self, other: T) -> Self {
        self.fmap(|x| x + other.into())
    }
}

/// Element-wise addition between &ComplexMatrix & Complex<f64>
impl<'a, T> Add<T> for &'a ComplexMatrix
where
    T: Into<Complex<f64>> + Copy,
{
    type Output = ComplexMatrix;

    fn add(self, other: T) -> Self::Output {
        self.fmap(|x| x + other.into())
    }
}

// Element-wise addition between Complex<f64> & ComplexMatrix
///
/// # Examples
///
/// ```rust
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
///
/// fn main() {
///     let mut a = ml_complex_matrix("1.0+1.0i 2.0+2.0i;
///                                    4.0+4.0i 5.0+5.0i");
///     let a_exp = ml_complex_matrix("2.0+2.0i 3.0+3.0i;
///                                    5.0+5.0i 6.0+6.0i");
///     assert_eq!(Complex64::new(1_f64, 1_f64) + a, a_exp);
/// }
/// ```
impl Add<ComplexMatrix> for Complex<f64> {
    type Output = ComplexMatrix;

    fn add(self, other: ComplexMatrix) -> Self::Output {
        other.add(self)
    }
}

/// Element-wise addition between Complex<f64> & &ComplexMatrix
impl<'a> Add<&'a ComplexMatrix> for Complex<f64> {
    type Output = ComplexMatrix;

    fn add(self, other: &'a ComplexMatrix) -> Self::Output {
        other.add(self)
    }
}

// =============================================================================
// Standard Operation for Matrix (Neg)
// =============================================================================
/// Negation of Complex Matrix
///
/// # Examples
/// ```rust
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
///
/// let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
///                             Complex64::new(2f64, 2f64),
///                             Complex64::new(3f64, 3f64),
///                             Complex64::new(4f64, 4f64)],
///                             2, 2, Row);
/// let a_neg = complex_matrix(vec![Complex64::new(-1f64, -1f64),
///                                 Complex64::new(-2f64, -2f64),
///                                 Complex64::new(-3f64, -3f64),
///                                 Complex64::new(-4f64, -4f64)],
///                             2, 2, Row);
/// assert_eq!(-a, a_neg);
/// ```
impl Neg for ComplexMatrix {
    type Output = Self;

    fn neg(self) -> Self {
        complex_matrix(
            self.data
                .into_iter()
                .map(|x: Complex<f64>| -x)
                .collect::<Vec<Complex<f64>>>(),
            self.row,
            self.col,
            self.shape,
        )
    }
}

/// Negation of &'a Complex Matrix
impl<'a> Neg for &'a ComplexMatrix {
    type Output = ComplexMatrix;

    fn neg(self) -> Self::Output {
        complex_matrix(
            self.data
                .clone()
                .into_iter()
                .map(|x: Complex<f64>| -x)
                .collect::<Vec<Complex<f64>>>(),
            self.row,
            self.col,
            self.shape,
        )
    }
}

// =============================================================================
// Standard Operation for Matrix (Sub)
// =============================================================================
/// Subtraction between Complex Matrix
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
///     let a = ml_complex_matrix("10.0+10.0i 20.0+20.0i;
///                                40.0+40.0i 50.0+50.0i");
///     let b = ml_complex_matrix("1.0+1.0i 2.0+2.0i;
///                                4.0+4.0i 5.0+5.0i");
///     let diff = ml_complex_matrix("9.0+9.0i 18.0+18.0i;
///                                   36.0+36.0i 45.0+45.0i");
///     assert_eq!(a-b, diff);
/// }
/// ```
impl Sub<ComplexMatrix> for ComplexMatrix {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        assert_eq!(&self.row, &other.row);
        assert_eq!(&self.col, &other.col);
        let mut result = complex_matrix(self.data.clone(), self.row, self.col, self.shape);
        for i in 0..self.row {
            for j in 0..self.col {
                result[(i, j)] -= other[(i, j)];
            }
        }
        result
    }
}

impl<'a, 'b> Sub<&'b ComplexMatrix> for &'a ComplexMatrix {
    type Output = ComplexMatrix;

    fn sub(self, rhs: &'b ComplexMatrix) -> Self::Output {
        self.sub_vec(rhs)
    }
}

/// Subtraction between Complex Matrix & Complex<f64>
impl<T> Sub<T> for ComplexMatrix
where
    T: Into<Complex<f64>> + Copy,
{
    type Output = Self;

    fn sub(self, other: T) -> Self::Output {
        self.fmap(|x| x - other.into())
    }
}

/// Subtraction between &Complex Matrix & Complex<f64>
impl<'a, T> Sub<T> for &'a ComplexMatrix
where
    T: Into<Complex<f64>> + Copy,
{
    type Output = ComplexMatrix;

    fn sub(self, other: T) -> Self::Output {
        self.fmap(|x| x - other.into())
    }
}

/// Subtraction Complex Matrix with Complex<f64>
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
///     let mut a = ml_complex_matrix("1.0+1.0i 2.0+2.0i;
///                                    4.0+4.0i 5.0+5.0i");
///     let a_exp = ml_complex_matrix("0.0+0.0i 1.0+1.0i;
///                                    3.0+3.0i 4.0+4.0i");
///     assert_eq!(a - Complex64::new(1_f64, 1_f64), a_exp);
/// }
/// ```
impl Sub<ComplexMatrix> for Complex<f64> {
    type Output = ComplexMatrix;

    fn sub(self, other: ComplexMatrix) -> Self::Output {
        -other.sub(self)
    }
}

impl<'a> Sub<&'a ComplexMatrix> for f64 {
    type Output = ComplexMatrix;

    fn sub(self, other: &'a ComplexMatrix) -> Self::Output {
        -other.sub(self)
    }
}

// =============================================================================
// Multiplication for Complex Matrix
// =============================================================================
/// Element-wise multiplication between Complex Matrix vs Complex<f64>
impl Mul<Complex<f64>> for ComplexMatrix {
    type Output = Self;

    fn mul(self, other: Complex<f64>) -> Self::Output {
        self.fmap(|x| x * other)
    }
}

impl Mul<ComplexMatrix> for Complex<f64> {
    type Output = ComplexMatrix;

    fn mul(self, other: ComplexMatrix) -> Self::Output {
        other.mul(self)
    }
}

impl<'a> Mul<&'a ComplexMatrix> for Complex<f64> {
    type Output = ComplexMatrix;

    fn mul(self, other: &'a ComplexMatrix) -> Self::Output {
        other.mul_scalar(self)
    }
}

/// Matrix Multiplication
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
///     let mut a = ml_complex_matrix("1.0+1.0i 2.0+2.0i;
///                                    4.0+4.0i 5.0+5.0i");
///     let mut b = ml_complex_matrix("2.0+2.0i 2.0+2.0i;
///                                    5.0+5.0i 5.0+5.0i");
///     let prod = ml_complex_matrix("0.0+24.0i 0.0+24.0i;
///                                    0.0+66.0i 0.0+66.0i");
///     assert_eq!(a * b, prod);
/// }
/// ```
impl Mul<ComplexMatrix> for ComplexMatrix {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        matmul(&self, &other)
    }
}

impl<'a, 'b> Mul<&'b ComplexMatrix> for &'a ComplexMatrix {
    type Output = ComplexMatrix;

    fn mul(self, other: &'b ComplexMatrix) -> Self::Output {
        matmul(self, other)
    }
}

#[allow(non_snake_case)]
impl Mul<Vec<Complex<f64>>> for ComplexMatrix {
    type Output = Vec<Complex<f64>>;

    fn mul(self, other: Vec<Complex<f64>>) -> Self::Output {
        self.apply(&other)
    }
}

#[allow(non_snake_case)]
impl<'a, 'b> Mul<&'b Vec<Complex<f64>>> for &'a ComplexMatrix {
    type Output = Vec<Complex<f64>>;

    fn mul(self, other: &'b Vec<Complex<f64>>) -> Self::Output {
        self.apply(other)
    }
}

/// Matrix multiplication for `Vec<Complex<f64>>` vs `ComplexMatrix`
impl Mul<ComplexMatrix> for Vec<Complex<f64>> {
    type Output = Vec<Complex<f64>>;

    fn mul(self, other: ComplexMatrix) -> Self::Output {
        assert_eq!(self.len(), other.row);
        let mut c = vec![Complex::zero(); other.col];
        complex_gevm(Complex::one(), &self, &other, Complex::zero(), &mut c);
        c
    }
}

impl<'a, 'b> Mul<&'b ComplexMatrix> for &'a Vec<Complex<f64>> {
    type Output = Vec<Complex<f64>>;

    fn mul(self, other: &'b ComplexMatrix) -> Self::Output {
        assert_eq!(self.len(), other.row);
        let mut c = vec![Complex::zero(); other.col];
        complex_gevm(Complex::one(), self, other, Complex::zero(), &mut c);
        c
    }
}

// =============================================================================
// Standard Operation for Matrix (DIV)
// =============================================================================
/// Element-wise division between Complex Matrix vs Complex<f64>
impl Div<Complex<f64>> for ComplexMatrix {
    type Output = Self;

    fn div(self, other: Complex<f64>) -> Self::Output {
        self.fmap(|x| x / other)
    }
}

impl<'a> Div<Complex<f64>> for &'a ComplexMatrix {
    type Output = ComplexMatrix;

    fn div(self, other: Complex<f64>) -> Self::Output {
        self.fmap(|x| x / other)
    }
}

/// Index for Complex Matrix
///
/// `(usize, usize) -> Complex<f64>`
///
/// # Examples
/// ```rust
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
///
/// let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
///                             Complex64::new(2f64, 2f64),
///                             Complex64::new(3f64, 3f64),
///                             Complex64::new(4f64, 4f64)],
///                             2, 2, Row
///     );
/// assert_eq!(a[(0,1)], Complex64::new(2f64, 2f64));
/// ```
impl Index<(usize, usize)> for ComplexMatrix {
    type Output = Complex<f64>;

    fn index(&self, pair: (usize, usize)) -> &Complex<f64> {
        let p = self.ptr();
        let i = pair.0;
        let j = pair.1;
        assert!(i < self.row && j < self.col, "Index out of range");
        match self.shape {
            Shape::Row => unsafe { &*p.add(i * self.col + j) },
            Shape::Col => unsafe { &*p.add(i + j * self.row) },
        }
    }
}

/// IndexMut for Complex Matrix (Assign)
///
/// `(usize, usize) -> Complex<f64>`
///
/// # Examples
/// ```rust
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
///
/// let mut a = complex_matrix(vec![Complex64::new(1f64, 1f64),
///                             Complex64::new(2f64, 2f64),
///                             Complex64::new(3f64, 3f64),
///                             Complex64::new(4f64, 4f64)],
///                             2, 2, Row
///     );
/// assert_eq!(a[(0,1)], Complex64::new(2f64, 2f64));
/// ```
impl IndexMut<(usize, usize)> for ComplexMatrix {
    fn index_mut(&mut self, pair: (usize, usize)) -> &mut Complex<f64> {
        let i = pair.0;
        let j = pair.1;
        let r = self.row;
        let c = self.col;
        assert!(i < self.row && j < self.col, "Index out of range");
        let p = self.mut_ptr();
        match self.shape {
            Shape::Row => {
                let idx = i * c + j;
                unsafe { &mut *p.add(idx) }
            }
            Shape::Col => {
                let idx = i + j * r;
                unsafe { &mut *p.add(idx) }
            }
        }
    }
}

// =============================================================================
// Functional Programming Tools (Hand-written)
// =============================================================================

impl FPMatrix for ComplexMatrix {
    type Scalar = Complex<f64>;

    fn take_row(&self, n: usize) -> Self {
        if n >= self.row {
            return self.clone();
        }
        match self.shape {
            Shape::Row => {
                let new_data = self
                    .data
                    .clone()
                    .into_iter()
                    .take(n * self.col)
                    .collect::<Vec<Complex<f64>>>();
                complex_matrix(new_data, n, self.col, Shape::Row)
            }
            Shape::Col => {
                let mut temp_data: Vec<Complex<f64>> = Vec::new();
                for i in 0..n {
                    temp_data.extend(self.row(i));
                }
                complex_matrix(temp_data, n, self.col, Shape::Row)
            }
        }
    }

    fn take_col(&self, n: usize) -> Self {
        if n >= self.col {
            return self.clone();
        }
        match self.shape {
            Shape::Col => {
                let new_data = self
                    .data
                    .clone()
                    .into_iter()
                    .take(n * self.row)
                    .collect::<Vec<Complex<f64>>>();
                complex_matrix(new_data, self.row, n, Shape::Col)
            }
            Shape::Row => {
                let mut temp_data: Vec<Complex<f64>> = Vec::new();
                for i in 0..n {
                    temp_data.extend(self.col(i));
                }
                complex_matrix(temp_data, self.row, n, Shape::Col)
            }
        }
    }

    fn skip_row(&self, n: usize) -> Self {
        assert!(n < self.row, "Skip range is larger than row of matrix");

        let mut temp_data: Vec<Complex<f64>> = Vec::new();
        for i in n..self.row {
            temp_data.extend(self.row(i));
        }
        complex_matrix(temp_data, self.row - n, self.col, Shape::Row)
    }

    fn skip_col(&self, n: usize) -> Self {
        assert!(n < self.col, "Skip range is larger than col of matrix");

        let mut temp_data: Vec<Complex<f64>> = Vec::new();
        for i in n..self.col {
            temp_data.extend(self.col(i));
        }
        complex_matrix(temp_data, self.row, self.col - n, Shape::Col)
    }

    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Complex<f64>) -> Complex<f64>,
    {
        let result = self
            .data
            .iter()
            .map(|x| f(*x))
            .collect::<Vec<Complex<f64>>>();
        complex_matrix(result, self.row, self.col, self.shape)
    }

    /// Column map
    ///
    /// # Example
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    /// use peroxide::traits::fp::FPMatrix;
    ///
    /// fn main() {
    ///     let x = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    ///     let y = x.col_map(|r| r.fmap(|t| t + r[0]));
    ///
    ///     let y_col_map = complex_matrix(vec![Complex64::new(2f64, 2f64),
    ///                                         Complex64::new(4f64, 4f64),
    ///                                         Complex64::new(4f64, 4f64),
    ///                                         Complex64::new(6f64, 6f64)],
    ///                             2, 2, Col
    ///     );
    ///
    ///     assert_eq!(y, y_col_map);
    /// }
    /// ```
    fn col_map<F>(&self, f: F) -> ComplexMatrix
    where
        F: Fn(Vec<Complex<f64>>) -> Vec<Complex<f64>>,
    {
        let mut result = complex_matrix(
            vec![Complex::zero(); self.row * self.col],
            self.row,
            self.col,
            Shape::Col,
        );

        for i in 0..self.col {
            result.subs_col(i, &f(self.col(i)));
        }

        result
    }

    /// Row map
    ///
    /// # Example
    /// ```rust
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    /// use peroxide::traits::fp::FPMatrix;
    ///
    /// fn main() {
    ///     let x = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    ///     let y = x.row_map(|r| r.fmap(|t| t + r[0]));
    ///
    ///     let y_row_map = complex_matrix(vec![Complex64::new(2f64, 2f64),
    ///                                         Complex64::new(3f64, 3f64),
    ///                                         Complex64::new(6f64, 6f64),
    ///                                         Complex64::new(7f64, 7f64)],
    ///                             2, 2, Row
    ///     );
    ///
    ///     assert_eq!(y, y_row_map);
    /// }
    /// ```
    fn row_map<F>(&self, f: F) -> ComplexMatrix
    where
        F: Fn(Vec<Complex<f64>>) -> Vec<Complex<f64>>,
    {
        let mut result = complex_matrix(
            vec![Complex::zero(); self.row * self.col],
            self.row,
            self.col,
            Shape::Row,
        );

        for i in 0..self.row {
            result.subs_row(i, &f(self.row(i)));
        }

        result
    }

    fn col_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<Complex<f64>>) -> Vec<Complex<f64>>,
    {
        for i in 0..self.col {
            unsafe {
                let mut p = self.col_mut(i);
                let fv = f(self.col(i));
                for j in 0..p.len() {
                    *p[j] = fv[j];
                }
            }
        }
    }

    fn row_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<Complex<f64>>) -> Vec<Complex<f64>>,
    {
        for i in 0..self.col {
            unsafe {
                let mut p = self.row_mut(i);
                let fv = f(self.row(i));
                for j in 0..p.len() {
                    *p[j] = fv[j];
                }
            }
        }
    }

    fn reduce<F, T>(&self, init: T, f: F) -> Complex<f64>
    where
        F: Fn(Complex<f64>, Complex<f64>) -> Complex<f64>,
        T: Into<Complex<f64>>,
    {
        self.data.iter().fold(init.into(), |x, y| f(x, *y))
    }

    fn zip_with<F>(&self, f: F, other: &ComplexMatrix) -> Self
    where
        F: Fn(Complex<f64>, Complex<f64>) -> Complex<f64>,
    {
        assert_eq!(self.data.len(), other.data.len());
        let mut a = other.clone();
        if self.shape != other.shape {
            a = a.change_shape();
        }
        let result = self
            .data
            .iter()
            .zip(a.data.iter())
            .map(|(x, y)| f(*x, *y))
            .collect::<Vec<Complex<f64>>>();
        complex_matrix(result, self.row, self.col, self.shape)
    }

    fn col_reduce<F>(&self, f: F) -> Vec<Complex<f64>>
    where
        F: Fn(Vec<Complex<f64>>) -> Complex<f64>,
    {
        let mut v = vec![Complex::zero(); self.col];
        for i in 0..self.col {
            v[i] = f(self.col(i));
        }
        v
    }

    fn row_reduce<F>(&self, f: F) -> Vec<Complex<f64>>
    where
        F: Fn(Vec<Complex<f64>>) -> Complex<f64>,
    {
        let mut v = vec![Complex::zero(); self.row];
        for i in 0..self.row {
            v[i] = f(self.row(i));
        }
        v
    }
}

// =============================================================================
// Linear Algebra
// =============================================================================

/// Linear algebra trait
pub trait LinearAlgebra {
    fn back_subs(&self, b: &Vec<Complex<f64>>) -> Vec<Complex<f64>>;
    fn forward_subs(&self, b: &Vec<Complex<f64>>) -> Vec<Complex<f64>>;
    fn block(&self) -> (ComplexMatrix, ComplexMatrix, ComplexMatrix, ComplexMatrix);
    fn is_symmetric(&self) -> bool;
    // ToDo: Add other fn of this trait from src/structure/matrix.rs
}

impl LinearAlgebra for ComplexMatrix {
    /// Backward Substitution for Upper Triangular
    fn back_subs(&self, b: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let n = self.col;
        let mut y = vec![Complex::zero(); n];
        y[n - 1] = b[n - 1] / self[(n - 1, n - 1)];
        for i in (0..n - 1).rev() {
            let mut s = Complex::zero();
            for j in i + 1..n {
                s += self[(i, j)] * y[j];
            }
            y[i] = 1f64 / self[(i, i)] * (b[i] - s);
        }
        y
    }

    /// Forward substitution for Lower Triangular
    fn forward_subs(&self, b: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        let n = self.col;
        let mut y = vec![Complex::zero(); n];
        y[0] = b[0] / self[(0, 0)];
        for i in 1..n {
            let mut s = Complex::zero();
            for j in 0..i {
                s += self[(i, j)] * y[j];
            }
            y[i] = 1f64 / self[(i, i)] * (b[i] - s);
        }
        y
    }

    /// Block Partition
    ///
    /// # Examples
    /// ```rust
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    /// use num_complex::Complex64;
    /// use peroxide::complex::matrix::*;
    /// use peroxide::complex::matrix::LinearAlgebra;
    ///
    /// fn main() {
    ///     let a = complex_matrix(vec![Complex64::new(1f64, 1f64),
    ///                                 Complex64::new(2f64, 2f64),
    ///                                 Complex64::new(3f64, 3f64),
    ///                                 Complex64::new(4f64, 4f64)],
    ///                             2, 2, Row
    ///     );
    ///     let (m1, m2, m3, m4) = a.block();
    ///     assert_eq!(m1, ml_complex_matrix("1.0+1.0i"));
    ///     assert_eq!(m2, ml_complex_matrix("2.0+2.0i"));
    ///     assert_eq!(m3, ml_complex_matrix("3.0+3.0i"));
    ///     assert_eq!(m4, ml_complex_matrix("4.0+4.0i"));
    /// }
    /// ```
    fn block(&self) -> (Self, Self, Self, Self) {
        let r = self.row;
        let c = self.col;
        let l_r = self.row / 2;
        let l_c = self.col / 2;
        let r_l = r - l_r;
        let c_l = c - l_c;

        let mut m1 = complex_matrix(vec![Complex::zero(); l_r * l_c], l_r, l_c, self.shape);
        let mut m2 = complex_matrix(vec![Complex::zero(); l_r * c_l], l_r, c_l, self.shape);
        let mut m3 = complex_matrix(vec![Complex::zero(); r_l * l_c], r_l, l_c, self.shape);
        let mut m4 = complex_matrix(vec![Complex::zero(); r_l * c_l], r_l, c_l, self.shape);

        for idx_row in 0..r {
            for idx_col in 0..c {
                match (idx_row, idx_col) {
                    (i, j) if (i < l_r) && (j < l_c) => {
                        m1[(i, j)] = self[(i, j)];
                    }
                    (i, j) if (i < l_r) && (j >= l_c) => {
                        m2[(i, j - l_c)] = self[(i, j)];
                    }
                    (i, j) if (i >= l_r) && (j < l_c) => {
                        m3[(i - l_r, j)] = self[(i, j)];
                    }
                    (i, j) if (i >= l_r) && (j >= l_c) => {
                        m4[(i - l_r, j - l_c)] = self[(i, j)];
                    }
                    _ => (),
                }
            }
        }
        (m1, m2, m3, m4)
    }

    fn is_symmetric(&self) -> bool {
        if self.row != self.col {
            return false;
        }

        for i in 0..self.row {
            for j in i..self.col {
                if (!nearly_eq(self[(i, j)].re, self[(j, i)].re))
                    && (!nearly_eq(self[(i, j)].im, self[(j, i)].im))
                {
                    return false;
                }
            }
        }
        true
    }
}

impl MutMatrix for ComplexMatrix {
    type Scalar = Complex<f64>;

    unsafe fn col_mut(&mut self, idx: usize) -> Vec<*mut Complex<f64>> {
        assert!(idx < self.col, "Index out of range");
        match self.shape {
            Shape::Col => {
                let mut v: Vec<*mut Complex<f64>> = vec![&mut Complex::zero(); self.row];
                let start_idx = idx * self.row;
                let p = self.mut_ptr();
                for (i, j) in (start_idx..start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Row => {
                let mut v: Vec<*mut Complex<f64>> = vec![&mut Complex::zero(); self.row];
                let p = self.mut_ptr();
                for i in 0..v.len() {
                    v[i] = p.add(idx + i * self.col);
                }
                v
            }
        }
    }

    unsafe fn row_mut(&mut self, idx: usize) -> Vec<*mut Complex<f64>> {
        assert!(idx < self.row, "Index out of range");
        match self.shape {
            Shape::Row => {
                let mut v: Vec<*mut Complex<f64>> = vec![&mut Complex::zero(); self.col];
                let start_idx = idx * self.col;
                let p = self.mut_ptr();
                for (i, j) in (start_idx..start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Col => {
                let mut v: Vec<*mut Complex<f64>> = vec![&mut Complex::zero(); self.col];
                let p = self.mut_ptr();
                for i in 0..v.len() {
                    v[i] = p.add(idx + i * self.row);
                }
                v
            }
        }
    }

    unsafe fn swap(&mut self, idx1: usize, idx2: usize, shape: Shape) {
        match shape {
            Shape::Col => swap_complex_vec_ptr(&mut self.col_mut(idx1), &mut self.col_mut(idx2)),
            Shape::Row => swap_complex_vec_ptr(&mut self.row_mut(idx1), &mut self.row_mut(idx2)),
        }
    }

    unsafe fn swap_with_perm(&mut self, p: &Vec<(usize, usize)>, shape: Shape) {
        for (i, j) in p.iter() {
            self.swap(*i, *j, shape)
        }
    }
}

// ToDo: Move swap_complex_vec_ptr to low_level.rs
pub unsafe fn swap_complex_vec_ptr(
    lhs: &mut Vec<*mut Complex<f64>>,
    rhs: &mut Vec<*mut Complex<f64>>,
) {
    assert_eq!(lhs.len(), rhs.len(), "Should use same length vectors");
    for (&mut l, &mut r) in lhs.iter_mut().zip(rhs.iter_mut()) {
        std::ptr::swap(l, r);
    }
}

// =============================================================================
// Back-end Utils
// =============================================================================
/// Combine separated Complex Matrix to one Complex Matrix
///
/// # Examples
/// ```rust
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
/// use peroxide::traits::fp::FPMatrix;
///
/// fn main() {
///     let x1 = complex_matrix(vec![Complex64::new(1f64, 1f64)], 1, 1, Row);
///     let x2 = complex_matrix(vec![Complex64::new(2f64, 2f64)], 1, 1, Row);
///     let x3 = complex_matrix(vec![Complex64::new(3f64, 3f64)], 1, 1, Row);
///     let x4 = complex_matrix(vec![Complex64::new(4f64, 4f64)], 1, 1, Row);
///
///     let y = complex_combine(x1, x2, x3, x4);
///
///     let y_exp = complex_matrix(vec![Complex64::new(1f64, 1f64),
///                                     Complex64::new(2f64, 2f64),
///                                     Complex64::new(3f64, 3f64),
///                                     Complex64::new(4f64, 4f64)],
///                             2, 2, Row
///     );
///
///     assert_eq!(y, y_exp);
/// }
/// ```
pub fn complex_combine(
    m1: ComplexMatrix,
    m2: ComplexMatrix,
    m3: ComplexMatrix,
    m4: ComplexMatrix,
) -> ComplexMatrix {
    let l_r = m1.row;
    let l_c = m1.col;
    let c_l = m2.col;
    let r_l = m3.row;

    let r = l_r + r_l;
    let c = l_c + c_l;

    let mut m = complex_matrix(vec![Complex::zero(); r * c], r, c, m1.shape);

    for idx_row in 0..r {
        for idx_col in 0..c {
            match (idx_row, idx_col) {
                (i, j) if (i < l_r) && (j < l_c) => {
                    m[(i, j)] = m1[(i, j)];
                }
                (i, j) if (i < l_r) && (j >= l_c) => {
                    m[(i, j)] = m2[(i, j - l_c)];
                }
                (i, j) if (i >= l_r) && (j < l_c) => {
                    m[(i, j)] = m3[(i - l_r, j)];
                }
                (i, j) if (i >= l_r) && (j >= l_c) => {
                    m[(i, j)] = m4[(i - l_r, j - l_c)];
                }
                _ => (),
            }
        }
    }
    m
}

/// Inverse of Lower matrix
///
/// # Examples
///  ```rust
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
///
/// fn main() {
///     let a = ml_complex_matrix("2.0+2.0i 0.0+0.0i;
///                                2.0+2.0i 1.0+1.0i");
///     let b = complex_matrix(vec![Complex64::new(2f64, 2f64),
///                                 Complex64::new(0f64, 0f64),
///                                 Complex64::new(-2f64, -2f64),
///                                 Complex64::new(1f64, 1f64)],
///                             2, 2, Row
///     );
///     assert_eq!(complex_inv_l(a), b);
/// }
/// ```
pub fn complex_inv_l(l: ComplexMatrix) -> ComplexMatrix {
    let mut m = l.clone();

    match l.row {
        1 => l,
        2 => {
            m[(1, 0)] = -m[(1, 0)];
            m
        }
        _ => {
            let (l1, l2, l3, l4) = l.block();

            let m1 = complex_inv_l(l1);
            let m2 = l2;
            let m4 = complex_inv_l(l4);
            let m3 = -(&(&m4 * &l3) * &m1);

            complex_combine(m1, m2, m3, m4)
        }
    }
}

/// Inverse of upper triangular matrix
///
/// # Examples
///  ```rust
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
/// use peroxide::complex::matrix::*;
///
/// fn main() {
///     let a = ml_complex_matrix("2.0+2.0i 2.0+2.0i;
///                                0.0+0.0i 1.0+1.0i");
///     let b = complex_matrix(vec![Complex64::new(0.25f64, -0.25f64),
///                                 Complex64::new(-0.5f64, 0.5f64),
///                                 Complex64::new(0.0f64, 0.0f64),
///                                 Complex64::new(0.5f64, -0.5f64)],
///                             2, 2, Row
///     );
///     assert_eq!(complex_inv_u(a), b);
/// }
/// ```
pub fn complex_inv_u(u: ComplexMatrix) -> ComplexMatrix {
    let mut w = u.clone();

    match u.row {
        1 => {
            w[(0, 0)] = 1f64 / w[(0, 0)];
            w
        }
        2 => {
            let a = w[(0, 0)];
            let b = w[(0, 1)];
            let c = w[(1, 1)];
            let d = a * c;

            w[(0, 0)] = 1f64 / a;
            w[(0, 1)] = -b / d;
            w[(1, 1)] = 1f64 / c;
            w
        }
        _ => {
            let (u1, u2, u3, u4) = u.block();
            let m1 = complex_inv_u(u1);
            let m3 = u3;
            let m4 = complex_inv_u(u4);
            let m2 = -(m1.clone() * u2 * m4.clone());

            complex_combine(m1, m2, m3, m4)
        }
    }
}

/// Matrix multiply back-ends
pub fn matmul(a: &ComplexMatrix, b: &ComplexMatrix) -> ComplexMatrix {
    assert_eq!(a.col, b.row);
    let mut c = complex_matrix(vec![Complex::zero(); a.row * b.col], a.row, b.col, a.shape);
    complex_gemm(Complex::one(), a, b, Complex::zero(), &mut c);
    c
}

/// GEMM wrapper for Matrixmultiply
///
/// # Examples
/// ```rust
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use num_complex::Complex64;
///
/// use peroxide::complex::matrix::*;
///
/// fn main() {
///     let a = ml_complex_matrix("1.0+1.0i 2.0+2.0i;
///                                0.0+0.0i 1.0+1.0i");
///     let b = ml_complex_matrix("1.0+1.0i 0.0+0.0i;
///                                2.0+2.0i 1.0+1.0i");
///     let mut c1 = ml_complex_matrix("1.0+1.0i 1.0+1.0i;
///                                    1.0+1.0i 1.0+1.0i");
///     let mul_val = ml_complex_matrix("-10.0+10.0i -4.0+4.0i;
///                                      -4.0+4.0i -2.0+2.0i");
///
///     complex_gemm(Complex64::new(1.0, 1.0), &a, &b, Complex64::new(0.0, 0.0), &mut c1);
///     assert_eq!(c1, mul_val);
/// }
pub fn complex_gemm(
    alpha: Complex<f64>,
    a: &ComplexMatrix,
    b: &ComplexMatrix,
    beta: Complex<f64>,
    c: &mut ComplexMatrix,
) {
    let m = a.row;
    let k = a.col;
    let n = b.col;
    let (rsa, csa) = match a.shape {
        Shape::Row => (a.col as isize, 1isize),
        Shape::Col => (1isize, a.row as isize),
    };
    let (rsb, csb) = match b.shape {
        Shape::Row => (b.col as isize, 1isize),
        Shape::Col => (1isize, b.row as isize),
    };
    let (rsc, csc) = match c.shape {
        Shape::Row => (c.col as isize, 1isize),
        Shape::Col => (1isize, c.row as isize),
    };

    unsafe {
        matrixmultiply::zgemm(
            // Requires crate feature "cgemm"
            CGemmOption::Standard,
            CGemmOption::Standard,
            m,
            k,
            n,
            [alpha.re, alpha.im],
            a.ptr() as *const _,
            rsa,
            csa,
            b.ptr() as *const _,
            rsb,
            csb,
            [beta.re, beta.im],
            c.mut_ptr() as *mut _,
            rsc,
            csc,
        )
    }
}

/// General Matrix-Vector multiplication
pub fn complex_gemv(
    alpha: Complex<f64>,
    a: &ComplexMatrix,
    b: &Vec<Complex<f64>>,
    beta: Complex<f64>,
    c: &mut Vec<Complex<f64>>,
) {
    let m = a.row;
    let k = a.col;
    let n = 1usize;
    let (rsa, csa) = match a.shape {
        Shape::Row => (a.col as isize, 1isize),
        Shape::Col => (1isize, a.row as isize),
    };
    let (rsb, csb) = (1isize, 1isize);
    let (rsc, csc) = (1isize, 1isize);

    unsafe {
        matrixmultiply::zgemm(
            // Requires crate feature "cgemm"
            CGemmOption::Standard,
            CGemmOption::Standard,
            m,
            k,
            n,
            [alpha.re, alpha.im],
            a.ptr() as *const _,
            rsa,
            csa,
            b.as_ptr() as *const _,
            rsb,
            csb,
            [beta.re, beta.im],
            c.as_mut_ptr() as *mut _,
            rsc,
            csc,
        )
    }
}

/// General Vector-Matrix multiplication
pub fn complex_gevm(
    alpha: Complex<f64>,
    a: &Vec<Complex<f64>>,
    b: &ComplexMatrix,
    beta: Complex<f64>,
    c: &mut Vec<Complex<f64>>,
) {
    let m = 1usize;
    let k = a.len();
    let n = b.col;
    let (rsa, csa) = (1isize, 1isize);
    let (rsb, csb) = match b.shape {
        Shape::Row => (b.col as isize, 1isize),
        Shape::Col => (1isize, b.row as isize),
    };
    let (rsc, csc) = (1isize, 1isize);

    unsafe {
        matrixmultiply::zgemm(
            // Requires crate feature "cgemm"
            CGemmOption::Standard,
            CGemmOption::Standard,
            m,
            k,
            n,
            [alpha.re, alpha.im],
            a.as_ptr() as *const _,
            rsa,
            csa,
            b.ptr() as *const _,
            rsb,
            csb,
            [beta.re, beta.im],
            c.as_mut_ptr() as *mut _,
            rsc,
            csc,
        )
    }
}
