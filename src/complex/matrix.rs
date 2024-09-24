use std::{
    fmt,
    ops::{Add, Index, IndexMut, Neg, Sub},
};

use num_complex::Complex;
use rand_distr::num_traits::{zero, One, Zero};

use crate::{
    fuga::{
        gemv, InnerProduct, LinearOp, MatrixProduct, Normed, Series, Shape, TypedVector, Vector,
    },
    traits::fp::FPComplexMatrix,
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
            result[i] = self.row(i); // Sen: needs row implementation
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
        unimplemented!()
    }
}

impl Normed for ComplexMatrix {
    type UnsignedScalar = Complex<f64>;

    fn norm(&self, kind: crate::fuga::Norm) -> Self::UnsignedScalar {
        unimplemented!()
    }

    fn normalize(&self, kind: crate::fuga::Norm) -> Self
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
        // gemv(Complex::one(), self, other, Complex::zero(), &mut c);  // Todo: Sen
        c
    }
}

impl MatrixProduct for ComplexMatrix {
    fn kronecker(&self, other: &Self) -> Self {
        unimplemented!()
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
///     let mut a = ml_complex_matrix("10.0+10.0i 20.0+20.0i;
///                                    40.0+40.0i 50.0+50.0i");
///     let b = ml_complex_matrix("1.0+1.0i 2.0+2.0i;
///                                    4.0+4.0i 5.0+5.0i");
///     let diff = ml_complex_matrix("9.0+9.0i 18.0+18.0i;
///                                    36.0+36.0i 45.0+45.0i");
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
// Standard Operation for Matrix (DIV)
// =============================================================================

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

impl FPComplexMatrix for ComplexMatrix {
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
    /// use peroxide::traits::fp::FPComplexMatrix;
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
    /// use peroxide::traits::fp::FPComplexMatrix;
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
        unimplemented!()
    }

    fn row_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<Complex<f64>>) -> Vec<Complex<f64>>,
    {
        unimplemented!()
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
