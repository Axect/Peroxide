//! Matrix for Scientific computation
//!
//! ## Declare matrix
//!
//! * You can declare matrix by various ways.
//!     * R's way - Default
//!     * MATLAB's way
//!     * Python's way
//!     * Other macro
//!
//! ### R's way
//!
//! * Description: Same as R - `matrix(Vec<f64>, Row, Col, Shape)`
//! * Type: `matrix(Vec<T>, usize, usize, Shape) where T: std::convert::Into<f64> + Copy`
//!     * `Shape`: `Enum` for matrix shape - `Row` & `Col`
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix(c!(1,2,3,4), 2, 2, Row);
//!         a.print();
//!         //       c[0] c[1]
//!         // r[0]     1    2
//!         // r[1]     3    4
//!
//!         let b = matrix(c!(1,2,3,4), 2, 2, Col);
//!         b.print();
//!         //       c[0] c[1]
//!         // r[0]     1    3
//!         // r[1]     2    4
//!     }
//!     ```
//!
//! ### MATLAB's way
//!
//! * Description: Similar to MATLAB (But should use `&str`)
//! * Type: `ml_matrix(&str)`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = ml_matrix("1 2; 3 4");
//!         a.print();
//!         //       c[0] c[1]
//!         // r[0]     1    2
//!         // r[1]     3    4
//!     }
//!     ```
//!
//! ### Python's way
//!
//! * Description: Declare matrix as vector of vectors.
//! * Type: `py_matrix(Vec<Vec<T>>) where T: std::convert::Into<f64> + Copy`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = py_matrix(vec![vec![1, 2], vec![3, 4]]);
//!         a.print();
//!         //       c[0] c[1]
//!         // r[0]     1    2
//!         // r[1]     3    4
//!     }
//!     ```
//!
//! ### Other macro
//!
//! * Description: R-like macro to declare matrix
//! * For `R`,
//!
//!     ```R
//!     # R
//!     a = matrix(1:4, nrow = 2, ncol = 2, byrow = T)
//!     print(a)
//!     #      [,1] [,2]
//!     # [1,]    1    2
//!     # [2,]    3    4
//!     ```
//!
//! * For `Peroxide`,
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         a.print();
//!         //       c[0] c[1]
//!         // r[0]     1    2
//!         // r[1]     3    4
//!     }
//!     ```
//!
//! ## Basic Method for Matrix
//!
//! There are some useful methods for `Matrix`
//!
//! * `row(&self, index: usize) -> Vec<f64>` : Extract specific row as `Vec<f64>`
//! * `col(&self, index: usize) -> Vec<f64>` : Extract specific column as `Vec<f64>`
//! * `diag(&self) -> Vec<f64>`: Extract diagonal components as `Vec<f64>`
//! * `swap(&self, usize, usize, Shape)`: Swap two rows or columns (unsafe function)
//! * `subs_col(&mut self, usize, Vec<f64>)`: Substitute column with `Vec<f64>`
//! * `subs_row(&mut self, usize, Vec<f64>)`: Substitute row with `Vec<f64>`
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let mut a = ml_matrix("1 2; 3 4");
//!
//!         a.row(0).print(); // [1, 2]
//!         a.col(0).print(); // [1, 3]
//!         a.diag().print(); // [1, 4]
//!         unsafe {
//!             a.swap(0, 1, Row);
//!         }
//!         a.print();
//!         //      c[0] c[1]
//!         // r[0]    3    4
//!         // r[1]    1    2
//!
//!         let mut b = ml_matrix("1 2;3 4");
//!         b.subs_col(0, &c!(5, 6));
//!         b.subs_row(1, &c!(7, 8));
//!         b.print();
//!         //       c[0] c[1]
//!         // r[0]    5    2
//!         // r[1]    7    8
//!     }
//!     ```
//!
//! ## Read & Write
//!
//! In peroxide, we can write matrix to `csv`
//!
//! ### CSV (Not recommended)
//!
//! * `csv` feature should be required
//! * `write(&self, file_path: &str)`: Write matrix to csv
//! * `write_with_header(&self, file_path, header: Vec<&str>)`: Write with header
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!     # #[cfg(feature="csv")]
//!     # {
//!         let a = ml_matrix("1 2;3 4");
//!         a.write("example_data/matrix.csv").expect("Can't write file");
//!
//!         let b = ml_matrix("1 2; 3 4; 5 6");
//!         b.write_with_header("example_data/header.csv", vec!["odd", "even"])
//!             .expect("Can't write header file");
//!     # }
//!         println!("Complete!")
//!     }
//!     ```
//!
//! Also, you can read matrix from csv.
//!
//! * Type: `read(&str, bool, char) -> Result<Matrix, Box<Error>>`
//! * Description: `read(file_path, is_header, delimiter)`
//!
//!     ```no_run
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         # #[cfg(feature="csv")]
//!         # {
//!         // `csv` feature should be required
//!         let a = Matrix::read("example_data/matrix.csv", false, ',')
//!             .expect("Can't read matrix.csv file");
//!         a.print();
//!         # }
//!         //       c[0] c[1]
//!         // r[0]     1    2
//!         // r[1]     3    4
//!     }
//!     ```
//!
//! ### Convert to DataFrame (Recommended)
//!
//! To write columns or rows, `DataFrame` and `nc` feature could be the best choice.
//!
//! * `nc` feature should be required - `netcdf` or `libnetcdf` are prerequisites.
//! 
//! ```no_run
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = matrix(c!(1,2,3,4,5,6), 3, 2, Col);
//!
//!     // Construct DataFrame
//!     let mut df = DataFrame::new(vec![]);
//!     df.push("x", Series::new(a.col(0)));
//!     df.push("y", Series::new(a.col(1)));
//!
//!     // Write nc file (`nc` feature should be required)
//!     # #[cfg(feature="nc")]
//!     # {
//!     df.write_nc("data.nc").expect("Can't write data.nc");
//!
//!     // Read nc file (`nc` feature should be required)
//!     let dg = DataFrame::read_nc("data.nc").expect("Can't read data.nc");
//!     let x: Vec<f64> = dg["x"].to_vec();
//!     let y: Vec<f64> = dg["y"].to_vec();
//!
//!     assert_eq!(a.col(0), x);
//!     assert_eq!(a.col(1), y);
//!     # }
//! }
//! ```
//!
//! ## Concatenation
//!
//! There are two options to concatenate matrices.
//!
//! * `cbind`: Concatenate two matrices by column direction.
//! * `rbind`: Concatenate two matrices by row direction.
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = ml_matrix("1 2;3 4");
//!         let b = ml_matrix("5 6;7 8");
//!
//!         cbind(a.clone(), b.clone()).print();
//!         //      c[0] c[1] c[2] c[3]
//!         // r[0]    1    2    5    7
//!         // r[1]    3    4    6    8
//!
//!         rbind(a, b).print();
//!         //      c[0] c[1]
//!         // r[0]    1    2
//!         // r[1]    3    4
//!         // r[2]    5    6
//!         // r[3]    7    8
//!     }
//!     ```
//!
//! ## Matrix operations
//!
//! * In peroxide, can use basic operations between matrices. I'll show you by examples.
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         (a.clone() + 1).print(); // -, *, / are also available
//!         //      c[0] c[1]
//!         // r[0]    2    3
//!         // r[1]    4    5
//!
//!         let b = matrix!(5;8;1, 2, 2, Row);
//!         (a.clone() + b.clone()).print(); // - is also available
//!         //      c[0] c[1]
//!         // r[0]    6    8
//!         // r[1]   10   12
//!
//!         (a.clone() * b.clone()).print(); // Matrix multiplication
//!         //      c[0] c[1]
//!         // r[0]   19   22
//!         // r[1]   43   50
//!     }
//!     ```
//!
//! * `clone` is too annoying - We can use reference operations!
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!     fn main() {
//!         let a = ml_matrix("1 2;3 4");
//!         let b = ml_matrix("5 6;7 8");
//!
//!         (&a + 1).print();
//!         (&a + &b).print();
//!         (&a - &b).print();
//!         (&a * &b).print();
//!     }
//!     ```
//!
//! ## Extract & modify components
//!
//! * In peroxide, matrix data is saved as linear structure.
//! * But you can use two-dimensional index to extract or modify components.
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let mut a = matrix!(1;4;1, 2, 2, Row);
//!         a[(0,0)].print();   // 1
//!         a[(0,0)] = 2f64;    // Modify component
//!         a.print();
//!         //       c[0] c[1]
//!         //  r[0]    2    2
//!         //  r[1]    3    4
//!     }
//!     ```
//!
//! ## Conversion to `Vec<f64>`
//!
//! * Just use `row` or `col` method (I already showed at Basic method section).
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         a.row(0).print(); // [1, 2]
//!         assert_eq!(a.row(0), vec![1f64, 2f64]);
//!     }
//!     ```
//!
//! ## Useful constructor
//!
//! * `zeros(usize, usize)`: Construct matrix which elements are all zero
//! * `eye(usize)`: Identity matrix
//! * `rand(usize, usize)`: Construct random uniform matrix (from 0 to 1)
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = zeros(2, 2);
//!         assert_eq!(a, ml_matrix("0 0;0 0"));
//!
//!         let b = eye(2);
//!         assert_eq!(b, ml_matrix("1 0;0 1"));
//!
//!         let c = rand(2, 2);
//!         c.print(); // Random 2x2 matrix
//!     }
//!     ```
//! # Linear Algebra
//!
//! ## Transpose
//!
//! * Caution: Transpose does not consume the original value.
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         a.transpose().print();
//!         // Or you can use shorter one
//!         a.t().print();
//!         //      c[0] c[1]
//!         // r[0]    1    3
//!         // r[1]    2    4
//!     }
//!     ```
//!
//! ## LU Decomposition
//!
//! * Peroxide uses **complete pivoting** for LU decomposition - Very stable
//! * Since there are lots of causes to generate error, you should use `Option`
//! * `lu` returns `PQLU`
//!     * `PQLU` has four field - `p`, `q`, `l` , `u`
//!     * `p` means row permutations
//!     * `q` means column permutations
//!     * `l` means lower triangular matrix
//!     * `u` menas upper triangular matrix
//! * The structure of `PQLU` is as follows:
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     #[derive(Debug, Clone)]
//!     pub struct PQLU {
//!         pub p: Vec<usize>,
//!         pub q: Vec<usize>,
//!         pub l: Matrix,
//!         pub u: Matrix,
//!     }
//!     ```
//!
//! * Example of LU decomposition:
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix(c!(1,2,3,4), 2, 2, Row);
//!         let pqlu = a.lu();
//!         let (p,q,l,u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);
//!         assert_eq!(p, vec![1]); // swap 0 & 1 (Row)
//!         assert_eq!(q, vec![1]); // swap 0 & 1 (Col)
//!         assert_eq!(l, matrix(c!(1,0,0.5,1),2,2,Row));
//!         //      c[0] c[1]
//!         // r[0]    1    0
//!         // r[1]  0.5    1
//!         assert_eq!(u, matrix(c!(4,3,0,-0.5),2,2,Row));
//!         //      c[0] c[1]
//!         // r[0]    4    3
//!         // r[1]    0 -0.5
//!     }
//!     ```
//!
//! ## Determinant
//!
//! * Peroxide uses LU decomposition to obtain determinant ($ \mathcal{O}(n^3) $)
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         assert_eq!(a.det(), -2f64);
//!     }
//!     ```
//!
//! ## Inverse matrix
//!
//! * Peroxide uses LU decomposition (via GECP) to obtain inverse matrix.
//! * It needs two sub functions - `inv_l`, `inv_u`
//!     * For inverse of `L, U`, I use block partitioning. For example, for lower triangular matrix :
//!     $$ \begin{aligned} L &= \begin{pmatrix} L_1 & \mathbf{0} \\\ L_2 & L_3 \end{pmatrix} \\\ L^{-1} &= \begin{pmatrix} L_1^{-1} & \mathbf{0} \\\ -L_3^{-1}L_2 L_1^{-1} & L_3^{-1} \end{pmatrix} \end{aligned} $$
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         a.inv().print();
//!         //      c[0] c[1]
//!         // r[0]   -2    1
//!         // r[1]  1.5 -0.5
//!     }
//!     ```
//!
//! ## Tips for LU, Det, Inverse
//!
//! * If you save `self.lu()` rather than the direct use of `self.det()` or `self.lu()` then you
//! can get better performance (via memoization)
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = ml_matrix("1 2;3 4");
//!     let pqlu = a.lu();  // Memoization of LU
//!     pqlu.det().print(); // Same as a.det() but do not need an additional LU
//!     pqlu.inv().print(); // Same as a.inv() but do not need an additional LU
//! }
//! ```
//!
//! ## QR Decomposition (`O3` feature only)
//!
//! * Use `dgeqrf` of LAPACK
//! * Return `QR` structure.
//!     
//!     ```ignore
//!     pub struct QR {
//!         pub q: Matrix,
//!         pub r: Matrix,
//!     }
//!     ```
//!
//! * Example
//!     
//!     ```ignore
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     let a = ml_matrix("
//!         1 -1 4;
//!         1 4 -2;
//!         1 4 2;
//!         1 -1 0
//!     ");
//!
//!     // QR decomposition
//!     let qr = a.qr();
//!
//!     qr.q().print();
//!     //         c[0]    c[1]    c[2]
//!     // r[0]    -0.5     0.5    -0.5
//!     // r[1]    -0.5 -0.5000  0.5000
//!     // r[2]    -0.5 -0.5000    -0.5
//!     // r[3]    -0.5     0.5  0.5000
//!
//!     qr.r().print();
//!     //      c[0] c[1] c[2]
//!     // r[0]   -2   -3   -2
//!     // r[1]    0   -5    2
//!     // r[2]    0    0   -4
//!     ```
//!
//! ## Singular Value Decomposition (`O3` feature only)
//!
//! * Use `dgesvd` of LAPACK
//! * Return `SVD` structure
//!
//!     ```no_run
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     pub struct SVD {
//!         pub s: Vec<f64>,
//!         pub u: Matrix,
//!         pub vt: Matrix,
//!     }
//!     ```
//!
//! * Example
//!
//!     ```
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = ml_matrix("3 2 2;2 3 -2");
//!         #[cfg(feature="O3")]
//!         {
//!             // Full SVD
//!             let svd = a.svd();
//!             assert!(eq_vec(&vec![5f64, 3f64], &svd.s, 1e-7));
//!
//!             // Or Truncated SVD
//!             let svd2 = svd.truncated();
//!             assert!(eq_vec(&vec![5f64, 3f64], &svd2.s, 1e-7));
//!         }
//!         a.print();
//!     }
//!     ```
//!
//! ## Cholesky Decomposition
//!
//! * Use `dpotrf` of LAPACK
//! * Return Matrix (But there can be panic! - Not symmetric or Not positive definite)
//! * Example
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = ml_matrix("1 2;2 5");
//!         #[cfg(feature = "O3")]
//!         {
//!             let u = a.cholesky(Upper);
//!             assert_eq!(u, ml_matrix("1 2;0 1"));
//!
//!             let l = a.cholesky(Lower);
//!             assert_eq!(l, ml_matrix("1 0;2 1"));
//!         }
//!         a.print();
//!     }
//!     ```
//!
//! ## Moore-Penrose Pseudo Inverse
//!
//! * $ X^\dagger = \left(X^T X\right)^{-1} X^T $
//! * For `O3` feature, peroxide use SVD to obtain pseudo inverse
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         let pinv_a = a.pseudo_inv();
//!         let inv_a = a.inv();
//!
//!         assert_eq!(inv_a, pinv_a); // Nearly equal (not actually equal)
//!     }
//!     ```

#[cfg(feature="csv")]
extern crate csv;

#[cfg(feature = "O3")]
extern crate blas;
#[cfg(feature = "O3")]
extern crate lapack;

#[cfg(feature="csv")]
use self::csv::{ReaderBuilder, StringRecord, WriterBuilder};
use ::matrixmultiply;
#[cfg(feature = "O3")]
use blas::{daxpy, dgemm, dgemv};
#[cfg(feature = "O3")]
use lapack::{dgecon, dgeqrf, dgetrf, dgetri, dgetrs, dorgqr, dgesvd, dpotrf};
#[cfg(feature = "O3")]
use std::f64::NAN;

pub use self::Shape::{Col, Row};
use crate::numerical::eigen::{eigen, EigenMethod};
use crate::traits::{
    general::Algorithm,
    fp::{FPMatrix, FPVector},
    math::{InnerProduct, LinearOp, MatrixProduct, Norm, Normed, Vector},
    mutable::MutMatrix,
};
use crate::util::{
    low_level::{swap_vec_ptr, copy_vec_ptr},
    non_macro::{cbind, eye, rbind, zeros},
    useful::{nearly_eq, tab},
};
use crate::structure::dataframe::{Series, TypedVector};
use std::cmp::{max, min};
use std::convert;
pub use std::error::Error;
use std::fmt;
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};
use crate::traits::sugar::ScalableMut;

pub type Perms = Vec<(usize, usize)>;

/// To select matrices' binding.
///
/// Row - Row binding
/// Col - Column binding
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = matrix(vec![1,2,3,4], 2, 2, Row);
/// let b = matrix(vec![1,2,3,4], 2, 2, Col);
/// println!("{}", a); // Similar to [[1,2],[3,4]]
/// println!("{}", b); // Similar to [[1,3],[2,4]]
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Shape {
    Col,
    Row,
}

/// Print for Shape
impl fmt::Display for Shape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let to_display = match self {
            Row => "Row",
            Col => "Col",
        };
        write!(f, "{}", to_display)
    }
}

impl Default for Shape {
    fn default() -> Self {
        Shape::Col
    }
}

/// R-like matrix structure
///
/// # Examples
///
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = Matrix {
///     data: vec![1f64,2f64,3f64,4f64],
///     row: 2,
///     col: 2,
///     shape: Row,
/// }; // [[1,2],[3,4]]
/// ```
#[derive(Debug, Clone, Default)]
pub struct Matrix {
    pub data: Vec<f64>,
    pub row: usize,
    pub col: usize,
    pub shape: Shape,
}

// =============================================================================
// Various matrix constructor
// =============================================================================

/// R-like matrix constructor
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix(c!(1,2,3,4), 2, 2, Row);
///     a.print(); // Print matrix
/// }
/// ```
pub fn matrix<T>(v: Vec<T>, r: usize, c: usize, shape: Shape) -> Matrix
where
    T: convert::Into<f64>,
{
    Matrix {
        data: v.into_iter().map(|t| t.into()).collect::<Vec<f64>>(),
        row: r,
        col: c,
        shape,
    }
}

/// R-like matrix constructor (Explicit ver.)
pub fn r_matrix<T>(v: Vec<T>, r: usize, c: usize, shape: Shape) -> Matrix
where
    T: convert::Into<f64>,
{
    matrix(v, r, c, shape)
}

/// Python-like matrix constructor
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = py_matrix(vec![c!(1,2), c!(3,4)]);
///     let b = matrix(c!(1,2,3,4), 2, 2, Row);
///     assert_eq!(a, b);
/// }
/// ```
pub fn py_matrix<T>(v: Vec<Vec<T>>) -> Matrix
where
    T: convert::Into<f64> + Copy,
{
    let r = v.len();
    let c = v[0].len();
    let mut data = vec![0f64; r * c];
    for i in 0..r {
        for j in 0..c {
            let idx = i * c + j;
            data[idx] = v[i][j].into();
        }
    }
    matrix(data, r, c, Row)
}

/// Matlab-like matrix constructor
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = ml_matrix("1 2; 3 4");
///     let b = matrix(c!(1,2,3,4), 2, 2, Row);
///     assert_eq!(a, b);
/// }
/// ```
pub fn ml_matrix(s: &str) -> Matrix where {
    let str_rows: Vec<&str> = s.split(';').collect();
    let r = str_rows.len();
    let str_data = str_rows
        .into_iter()
        .map(|x| x.trim().split(' ').collect::<Vec<&str>>())
        .collect::<Vec<Vec<&str>>>();
    let c = str_data[0].len();
    let data = str_data
        .into_iter()
        .flat_map(|t| {
            t.into_iter()
                .map(|x| x.parse::<f64>().unwrap())
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<f64>>();
    matrix(data, r, c, Row)
}

/// Pretty Print
impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

/// PartialEq implements
impl PartialEq for Matrix {
    fn eq(&self, other: &Matrix) -> bool {
        if self.shape == other.shape {
            self.data
                .clone()
                .into_iter()
                .zip(other.data.clone())
                .all(|(x, y)| nearly_eq(x, y))
                && self.row == other.row
        } else {
            self.eq(&other.change_shape())
        }
    }
}

/// Main matrix structure
#[allow(dead_code)]
impl Matrix {
    /// Raw pointer for `self.data`
    pub fn ptr(&self) -> *const f64 {
        &self.data[0] as *const f64
    }

    /// Raw mutable pointer for `self.data`
    pub fn mut_ptr(&mut self) -> *mut f64 {
        &mut self.data[0] as *mut f64
    }

    /// Slice of `self.data`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// let a = matrix(vec![1,2,3,4], 2, 2, Col);
    /// let b = a.as_slice();
    /// assert_eq!(b, &[1f64,2f64,3f64,4f64]);
    /// ```
    pub fn as_slice<'a>(&'a self) -> &'a [f64] {
        &self.data[..]
    }

    /// Mutable slice of `self.data`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// let mut a = matrix(vec![1,2,3,4], 2, 2, Col);
    /// let mut b = a.as_mut_slice();
    /// b[0] = 5f64;
    /// assert_eq!(b, &[5f64, 2f64, 3f64, 4f64]);
    /// assert_eq!(a, matrix(vec![5,2,3,4], 2, 2, Col));
    /// ```
    pub fn as_mut_slice<'a>(&'a mut self) -> &'a mut [f64] {
        &mut self.data[..]
    }

    /// Change Bindings
    ///
    /// `Row` -> `Col` or `Col` -> `Row`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// let a = matrix(vec![1,2,3,4],2,2,Row);
    /// assert_eq!(a.shape, Row);
    /// let b = a.change_shape();
    /// assert_eq!(b.shape, Col);
    /// ```
    pub fn change_shape(&self) -> Self {
        let r = self.row;
        let c = self.col;
        assert_eq!(r * c, self.data.len());
        let l = r * c - 1;
        let mut data: Vec<f64> = self.data.clone();
        let ref_data = &self.data;

        match self.shape {
            Row => {
                for i in 0..l {
                    let s = (i * c) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                matrix(data, r, c, Col)
            }
            Col => {
                for i in 0..l {
                    let s = (i * r) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                matrix(data, r, c, Row)
            }
        }
    }

    /// Change Bindings Mutably
    ///
    /// `Row` -> `Col` or `Col` -> `Row`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// let a = matrix(vec![1,2,3,4],2,2,Row);
    /// assert_eq!(a.shape, Row);
    /// let b = a.change_shape();
    /// assert_eq!(b.shape, Col);
    /// ```
    pub fn change_shape_mut(&mut self) {
        let r = self.row;
        let c = self.col;
        assert_eq!(r * c, self.data.len());
        let l = r * c - 1;
        let ref_data = self.data.clone();

        match self.shape {
            Row => {
                for i in 0..l {
                    let s = (i * c) % l;
                    self.data[i] = ref_data[s];
                }
                self.data[l] = ref_data[l];
                self.shape = Col;
            }
            Col => {
                for i in 0..l {
                    let s = (i * r) % l;
                    self.data[i] = ref_data[s];
                }
                self.data[l] = ref_data[l];
                self.shape = Row;
            }
        }
    }

    /// Spread data(1D vector) to 2D formatted String
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// let a = matrix(vec![1,2,3,4],2,2,Row);
    /// println!("{}", a.spread()); // same as println!("{}", a);
    /// // Result:
    /// //       c[0] c[1]
    /// // r[0]     1    3
    /// // r[1]     2    4
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
                "Result is too large to print - {}x{}\nonly print {}x{} parts:\n{}",
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
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix(c!(1,2,3,4), 2, 2, Row);
    ///     assert_eq!(a.col(0), c!(1,3));
    /// }
    /// ```
    pub fn col(&self, index: usize) -> Vec<f64> {
        assert!(index < self.col);
        let mut container: Vec<f64> = vec![0f64; self.row];
        for i in 0..self.row {
            container[i] = self[(i, index)];
        }
        container
    }

    /// Extract Row
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix(c!(1,2,3,4), 2, 2, Row);
    ///     assert_eq!(a.row(0), c!(1,2));
    /// }
    /// ```
    pub fn row(&self, index: usize) -> Vec<f64> {
        assert!(index < self.row);
        let mut container: Vec<f64> = vec![0f64; self.col];
        for i in 0..self.col {
            container[i] = self[(index, i)];
        }
        container
    }

    /// Extract diagonal components
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix!(1;4;1, 2, 2, Row);
    ///     assert_eq!(a.diag(), c!(1,4));
    /// }
    /// ```
    pub fn diag(&self) -> Vec<f64> {
        let mut container = vec![0f64; self.row];
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
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// let a = matrix(vec![1,2,3,4], 2, 2, Row);
    /// println!("{}", a); // [[1,3],[2,4]]
    /// ```
    pub fn transpose(&self) -> Self {
        match self.shape {
            Row => matrix(self.data.clone(), self.col, self.row, Col),
            Col => matrix(self.data.clone(), self.col, self.row, Row),
        }
    }

    /// R-like transpose function
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix!(1;4;1, 2, 2, Row);
    ///     assert_eq!(a.transpose(), a.t());
    /// }
    /// ```
    pub fn t(&self) -> Self {
        self.transpose()
    }

    /// Write to CSV
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    ///     a.write("example_data/test.csv");
    /// }
    /// ```
    #[cfg(feature="csv")]
    pub fn write(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r = self.row;
        let c = self.col;
        for i in 0..r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for j in 0..c {
                record[j] = self[(i, j)].to_string();
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    /// Write to CSV (with round option)
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    ///     a.write_round("example_data/test.csv", 0);
    /// }
    /// ```
    #[cfg(feature="csv")]
    pub fn write_round(&self, file_path: &str, round: usize) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r = self.row;
        let c = self.col;
        for i in 0..r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for j in 0..c {
                record[j] = format!("{:.*}", round, self[(i, j)]);
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    #[cfg(feature="csv")]
    pub fn write_with_header(
        &self,
        file_path: &str,
        header: Vec<&str>,
    ) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r = self.row;
        let c = self.col;
        assert_eq!(c, header.len());
        wtr.write_record(header)?;
        for i in 0..r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for j in 0..c {
                record[j] = self[(i, j)].to_string();
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    #[cfg(feature="csv")]
    pub fn write_with_header_round(
        &self,
        file_path: &str,
        header: Vec<&str>,
        round: usize,
    ) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r = self.row;
        let c = self.col;
        assert_eq!(c, header.len());
        wtr.write_record(header)?;
        for i in 0..r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for j in 0..c {
                record[j] = format!("{:.*}", round, self[(i, j)]);
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    /// Read from CSV
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    /// use std::process;
    ///
    /// fn main() {
    ///     let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    ///     a.write_round("example_data/test.csv", 0);
    ///
    ///     let b = Matrix::read("example_data/test.csv", false, ','); // header = false, delimiter = ','
    ///     match b {
    ///         Ok(mat) => println!("{}", mat),
    ///         Err(err) => {
    ///             println!("{}", err);
    ///             process::exit(1);
    ///         }
    ///     }
    /// }
    /// ```
    #[cfg(feature="csv")]
    pub fn read(file_path: &str, header: bool, delimiter: char) -> Result<Matrix, Box<dyn Error>> {
        let mut rdr = ReaderBuilder::new()
            .has_headers(header)
            .delimiter(delimiter as u8)
            .from_path(file_path)?;

        let records = rdr
            .records()
            .collect::<Result<Vec<StringRecord>, csv::Error>>()?;
        let result = records;

        let l_row = result.len();
        let l_col = result[0].len();

        let mut m = matrix(vec![0f64; l_row * l_col], l_row, l_col, Row);

        for i in 0..l_row {
            for j in 0..l_col {
                m[(i, j)] = result[i][j].parse::<f64>().unwrap();
            }
        }

        Ok(m)
    }

    /// Should check shape
    pub fn subs(&mut self, idx: usize, v: &Vec<f64>) {
        let p = &mut self.mut_ptr();
        match self.shape {
            Row => {
                let c = self.col;
                unsafe {
                    p.add(idx * c).copy_from(v.as_ptr(), c);
                }
            }
            Col => {
                let r = self.row;
                unsafe {
                    p.add(idx * r).copy_from(v.as_ptr(), r);
                }
            }
        }
    }

    /// Substitute Col
    #[inline]
    pub fn subs_col(&mut self, idx: usize, v: &Vec<f64>) {
        for i in 0..self.row {
            self[(i, idx)] = v[i];
        }
    }

    /// Substitute Row
    #[inline]
    pub fn subs_row(&mut self, idx: usize, v: &Vec<f64>) {
        for j in 0..self.col {
            self[(idx, j)] = v[j];
        }
    }

    /// From index operations
    pub fn from_index<F, G>(f: F, size: (usize, usize)) -> Matrix
    where
        F: Fn(usize, usize) -> G + Copy,
        G: Into<f64>,
    {
        let row = size.0;
        let col = size.1;

        let mut mat = matrix(vec![0f64; row * col], row, col, Row);

        for i in 0..row {
            for j in 0..col {
                mat[(i, j)] = f(i, j).into();
            }
        }
        mat
    }

    /// Matrix to `Vec<Vec<f64>>`
    ///
    /// To send `Matrix` to `inline-python`
    pub fn to_vec(&self) -> Vec<Vec<f64>> {
        let mut result = vec![vec![0f64; self.col]; self.row];
        for i in 0..self.row {
            result[i] = self.row(i);
        }
        result
    }

    pub fn to_diag(&self) -> Matrix {
        assert!(self.row == self.col, "Should be square matrix");
        let mut result = matrix(vec![0f64; self.row * self.col], self.row, self.col, Row);
        let diag = self.diag();
        for i in 0..self.row {
            result[(i, i)] = diag[i];
        }
        result
    }

    /// Submatrix
    ///
    /// # Description
    /// Return below elements of matrix to new matrix
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
    /// 
    /// fn main() {
    ///     let a = ml_matrix("1 2 3;4 5 6;7 8 9");
    ///     let b = ml_matrix("5 6;8 9");
    ///     let c = a.submat((1, 1), (2, 2));
    ///     assert_eq!(b, c);   
    /// }
    /// ```
    pub fn submat(&self, start: (usize, usize), end: (usize, usize)) -> Matrix {
        let row = end.0 - start.0 + 1;
        let col = end.1 - start.1 + 1;
        let mut result = matrix(vec![0f64; row * col], row, col, self.shape);
        for i in 0 .. row {
            for j in 0 .. col {
                result[(i, j)] = self[(start.0 + i, start.1 + j)];
            }
        }
        result
    }

    /// Substitute matrix to specific position
    ///
    /// # Description
    /// Substitute below elements of matrix
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
    /// 
    /// fn main() {
    ///     let mut a = ml_matrix("1 2 3;4 5 6;7 8 9");
    ///     let b = ml_matrix("1 2;3 4");
    ///     let c = ml_matrix("1 2 3;4 1 2;7 3 4");
    ///     a.subs_mat((1,1), (2,2), &b);
    ///     assert_eq!(a, c);       
    /// }
    /// ```
    pub fn subs_mat(&mut self, start: (usize, usize), end: (usize, usize), m: &Matrix) {
        let row = end.0 - start.0 + 1;
        let col = end.1 - start.1 + 1;
        for i in 0 .. row {
            for j in 0 .. col {
                self[(start.0 + i, start.1 + j)] = m[(i, j)];
            }
        }
    }

    /// Matrix from series
    ///
    /// # Example
    /// ```rust
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = Series::new(c!(1,2,3,4));
    ///     let m = Matrix::from_series(&a, 2, 2, Row);
    ///     
    ///     assert_eq!(m, matrix(c!(1,2,3,4), 2, 2, Row));
    /// }
    /// ```
    pub fn from_series(series: &Series, row: usize, col: usize, shape: Shape) -> Self {
        let v: Vec<f64> = series.to_vec();
        matrix(v, row, col, shape)
    }
}

// =============================================================================
// Mathematics for Matrix
// =============================================================================
impl Vector for Matrix {
    type Scalar = f64;

    fn add_vec(&self, other: &Self) -> Self {
        assert_eq!(self.row, other.row);
        assert_eq!(self.col, other.col);

        match () {
            #[cfg(feature = "O3")]
            () => {
                if self.shape != other.shape {
                    return self.add(&other.change_shape());
                }
                let x = &self.data;
                let mut y = other.data.clone();
                let n_i32 = x.len() as i32;
                let a_f64 = 1f64;
                unsafe {
                    daxpy(n_i32, a_f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => {
                let mut result = matrix(self.data.clone(), self.row, self.col, self.shape);
                for i in 0..self.row {
                    for j in 0..self.col {
                        result[(i, j)] += other[(i, j)];
                    }
                }
                result
            }
        }
    }

    fn sub_vec(&self, other: &Self) -> Self {
        assert_eq!(self.row, other.row);
        assert_eq!(self.col, other.col);
        match () {
            #[cfg(feature = "O3")]
            () => {
                if self.shape != other.shape {
                    return self.sub(&other.change_shape());
                }
                let x = &other.data;
                let mut y = self.data.clone();
                let n_i32 = x.len() as i32;
                let a_f64 = -1f64;
                unsafe {
                    daxpy(n_i32, a_f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => {
                let mut result = matrix(self.data.clone(), self.row, self.col, self.shape);
                for i in 0..self.row {
                    for j in 0..self.col {
                        result[(i, j)] -= other[(i, j)];
                    }
                }
                result
            }
        }
    }

    fn mul_scalar(&self, other: Self::Scalar) -> Self {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &self.data;
                let mut y = vec![0f64; x.len()];
                let a_f64 = other;
                let n_i32 = x.len() as i32;

                unsafe {
                    daxpy(n_i32, a_f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => {
                let scalar = other;
                self.fmap(|x| x * scalar)
            }
        }
    }
}

impl Normed for Matrix {
    type UnsignedScalar = f64;
    fn norm(&self, kind: Norm) -> f64 {
        match kind {
            Norm::F => {
                let mut s = 0f64;
                for i in 0..self.data.len() {
                    s += self.data[i].powi(2);
                }
                s.sqrt()
            }
            Norm::Lpq(p, q) => {
                let mut s = 0f64;
                for j in 0..self.col {
                    let mut s_row = 0f64;
                    for i in 0..self.row {
                        s_row += self[(i, j)].powi(p as i32);
                    }
                    s += s_row.powf(q as f64 / (p as f64));
                }
                s.powf(1f64 / (q as f64))
            }
            Norm::L1 => {
                let mut m = std::f64::MIN;
                match self.shape {
                    Row => self.change_shape().norm(Norm::L1),
                    Col => {
                        for c in 0..self.col {
                            let s = self.col(c).iter().sum();
                            if s > m {
                                m = s;
                            }
                        }
                        m
                    }
                }
            }
            Norm::LInf => {
                let mut m = std::f64::MIN;
                match self.shape {
                    Col => self.change_shape().norm(Norm::LInf),
                    Row => {
                        for r in 0..self.row {
                            let s = self.row(r).iter().sum();
                            if s > m {
                                m = s;
                            }
                        }
                        m
                    }
                }
            }
            Norm::L2 => {
                let a = &self.t() * self;
                let eig = eigen(&a, EigenMethod::Jacobi);
                eig.eigenvalue.norm(Norm::LInf)
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
impl InnerProduct for Matrix {
    fn dot(&self, rhs: &Self) -> f64 {
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
impl LinearOp<Vec<f64>, Vec<f64>> for Matrix {
    fn apply(&self, other: &Vec<f64>) -> Vec<f64> {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = other;
                let mut y = vec![0f64; self.row];
                let A = &self.data;
                let m_i32 = self.row as i32;
                let n_i32 = self.col as i32;
                match self.shape {
                    Row => unsafe {
                        dgemv(b'T', n_i32, m_i32, 1f64, A, n_i32, x, 1, 0f64, &mut y, 1);
                    },
                    Col => unsafe {
                        dgemv(b'N', m_i32, n_i32, 1f64, A, m_i32, x, 1, 0f64, &mut y, 1);
                    },
                }
                y
            }
            _ => {
                assert_eq!(self.col, other.len());
                let mut c = vec![0f64; self.row];
                gemv(1f64, self, other, 0f64, &mut c);
                c
            }
        }
    }
}

impl MatrixProduct for Matrix {
    fn kronecker(&self, other: &Self) -> Self {
        let r1 = self.row;
        let c1 = self.col;

        let mut result = self[(0, 0)] * other;

        for j in 1..c1 {
            let n = self[(0, j)] * other;
            result = cbind(result, n);
        }

        for i in 1..r1 {
            let mut m = self[(i, 0)] * other;
            for j in 1..c1 {
                let n = self[(i, j)] * other;
                m = cbind(m, n);
            }
            result = rbind(result, m);
        }
        result
    }

    fn hadamard(&self, other: &Self) -> Matrix {
        assert_eq!(self.row, other.row);
        assert_eq!(self.col, other.col);

        let r = self.row;
        let c = self.col;

        let mut m = matrix(vec![0f64; r * c], r, c, self.shape);
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
/// `Matrix` to `Vec<f64>`
impl Into<Vec<f64>> for Matrix {
    fn into(self) -> Vec<f64> {
        self.data
    }
}

/// `&Matrix` to `&Vec<f64>`
impl<'a> Into<&'a Vec<f64>> for &'a Matrix {
    fn into(self) -> &'a Vec<f64> {
        &self.data
    }
}

/// `Vec<f64>` to `Matrix`
impl Into<Matrix> for Vec<f64> {
    fn into(self) -> Matrix {
        let l = self.len();
        matrix(self, l, 1, Col)
    }
}

impl Into<Matrix> for &Vec<f64> {
    fn into(self) -> Matrix {
        let l = self.len();
        matrix(self.clone(), l, 1, Col)
    }
}

// =============================================================================
// Standard Operation for Matrix (ADD)
// =============================================================================

/// Element-wise addition of Matrix
///
/// # Caution
/// > You should remember ownership.
/// > If you use Matrix `a,b` then you can't use them after.
impl Add<Matrix> for Matrix {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        assert_eq!(&self.row, &other.row);
        assert_eq!(&self.col, &other.col);

        match () {
            #[cfg(feature = "O3")]
            () => {
                if self.shape != other.shape {
                    return self.add(other.change_shape());
                }
                let x = &self.data;
                let mut y = other.data.clone();
                let n_i32 = x.len() as i32;
                let a_f64 = 1f64;
                unsafe {
                    daxpy(n_i32, a_f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => {
                let mut result = matrix(self.data.clone(), self.row, self.col, self.shape);
                for i in 0..self.row {
                    for j in 0..self.col {
                        result[(i, j)] += other[(i, j)];
                    }
                }
                result
            }
        }
    }
}

impl<'a, 'b> Add<&'b Matrix> for &'a Matrix {
    type Output = Matrix;

    fn add(self, rhs: &'b Matrix) -> Self::Output {
        self.add_vec(rhs)
    }
}

/// Element-wise addition between Matrix & f64
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;4;1, 2, 2, Row);
///     assert_eq!(a + 1, matrix!(2;5;1, 2, 2, Row));
/// }
/// ```
impl<T> Add<T> for Matrix
where
    T: convert::Into<f64> + Copy,
{
    type Output = Self;
    fn add(self, other: T) -> Self {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &self.data;
                let mut y = vec![other.into(); x.len()];
                let n_i32 = x.len() as i32;
                unsafe {
                    daxpy(n_i32, 1f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => self.fmap(|x| x + other.into()),
        }
    }
}

/// Element-wise addition between &Matrix & f64
impl<'a, T> Add<T> for &'a Matrix
where
    T: convert::Into<f64> + Copy,
{
    type Output = Matrix;
    fn add(self, other: T) -> Self::Output {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &self.data;
                let mut y = vec![other.into(); x.len()];
                let n_i32 = x.len() as i32;
                unsafe {
                    daxpy(n_i32, 1f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => self.fmap(|x| x + other.into()),
        }
    }
}

/// Element-wise addition between f64 & matrix
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;4;1, 2, 2, Row);
///     assert_eq!(1f64 + a, matrix!(2;5;1, 2, 2, Row));
/// }
/// ```
impl Add<Matrix> for f64 {
    type Output = Matrix;

    fn add(self, other: Matrix) -> Matrix {
        other.add(self)
    }
}

/// Element-wise addition between f64 & &Matrix
impl<'a> Add<&'a Matrix> for f64 {
    type Output = Matrix;

    fn add(self, other: &'a Matrix) -> Self::Output {
        other.add(self)
    }
}

/// Element-wise addition between i32 & Matrix
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;4;1, 2, 2, Row);
///     assert_eq!(1 + a, matrix!(2;5;1, 2, 2, Row));
/// }
/// ```
impl Add<Matrix> for i32 {
    type Output = Matrix;

    fn add(self, other: Matrix) -> Matrix {
        other.add(self)
    }
}

/// Element-wise addition between i32 & &Matrix
impl<'a> Add<&'a Matrix> for i32 {
    type Output = Matrix;

    fn add(self, other: &'a Matrix) -> Self::Output {
        other.add(self)
    }
}

/// Element-wise addition between usize & matrix
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;4;1, 2, 2, Row);
///     assert_eq!(1 as usize + a, matrix!(2;5;1, 2, 2, Row));
/// }
/// ```
impl Add<Matrix> for usize {
    type Output = Matrix;

    fn add(self, other: Matrix) -> Matrix {
        other.add(self as f64)
    }
}

/// Element-wise addition between usize & &Matrix
impl<'a> Add<&'a Matrix> for usize {
    type Output = Matrix;

    fn add(self, other: &'a Matrix) -> Self::Output {
        other.add(self as f64)
    }
}

// =============================================================================
// Standard Operation for Matrix (Neg)
// =============================================================================
/// Negation of Matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = matrix(vec![1,2,3,4],2,2,Row);
/// println!("{}", -a); // [[-1,-2],[-3,-4]]
/// ```
impl Neg for Matrix {
    type Output = Self;

    fn neg(self) -> Self {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &self.data;
                let mut y = vec![0f64; x.len()];
                let n_i32 = x.len() as i32;
                unsafe {
                    daxpy(n_i32, -1f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => matrix(
                self.data.into_iter().map(|x: f64| -x).collect::<Vec<f64>>(),
                self.row,
                self.col,
                self.shape,
            ),
        }
    }
}

/// Negation of &'a Matrix
impl<'a> Neg for &'a Matrix {
    type Output = Matrix;

    fn neg(self) -> Self::Output {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &self.data;
                let mut y = vec![0f64; x.len()];
                let n_i32 = x.len() as i32;
                unsafe {
                    daxpy(n_i32, -1f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => matrix(
                self.data
                    .clone()
                    .into_iter()
                    .map(|x: f64| -x)
                    .collect::<Vec<f64>>(),
                self.row,
                self.col,
                self.shape,
            ),
        }
    }
}

// =============================================================================
// Standard Operation for Matrix (Sub)
// =============================================================================
/// Subtraction between Matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = matrix(vec![1,2,3,4],2,2,Row);
/// let b = matrix(vec![1,2,3,4],2,2,Col);
/// println!("{}", a - b); // [[0, -1], [1, 0]]
/// ```
impl Sub<Matrix> for Matrix {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        assert_eq!(&self.row, &other.row);
        assert_eq!(&self.col, &other.col);
        match () {
            #[cfg(feature = "O3")]
            () => {
                if self.shape != other.shape {
                    return self.sub(other.change_shape());
                }
                let x = &other.data;
                let mut y = self.data.clone();
                let n_i32 = x.len() as i32;
                let a_f64 = -1f64;
                unsafe {
                    daxpy(n_i32, a_f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => {
                let mut result = matrix(self.data.clone(), self.row, self.col, self.shape);
                for i in 0..self.row {
                    for j in 0..self.col {
                        result[(i, j)] -= other[(i, j)];
                    }
                }
                result
            }
        }
    }
}

impl<'a, 'b> Sub<&'b Matrix> for &'a Matrix {
    type Output = Matrix;

    fn sub(self, rhs: &'b Matrix) -> Matrix {
        self.sub_vec(rhs)
    }
}

/// Subtraction between Matrix & f64
impl<T> Sub<T> for Matrix
where
    T: convert::Into<f64> + Copy,
{
    type Output = Self;

    fn sub(self, other: T) -> Self {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let mut y = self.data;
                let x = vec![other.into(); y.len()];
                let n_i32 = y.len() as i32;
                unsafe {
                    daxpy(n_i32, -1f64, &x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => self.fmap(|x| x - other.into()),
        }
    }
}

/// Subtraction between &Matrix & f64
impl<'a, T> Sub<T> for &'a Matrix
where
    T: convert::Into<f64> + Copy,
{
    type Output = Matrix;

    fn sub(self, other: T) -> Self::Output {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let mut y = self.data.clone();
                let x = vec![other.into(); y.len()];
                let n_i32 = y.len() as i32;
                unsafe {
                    daxpy(n_i32, -1f64, &x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => self.fmap(|x| x - other.into()),
        }
    }
}

/// Subtraction Matrix with f64
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix(vec![1,2,3,4],2,2,Row);
///     assert_eq!(a - 1f64, matrix!(0;3;1, 2, 2, Row));
/// }
/// ```
impl Sub<Matrix> for f64 {
    type Output = Matrix;

    fn sub(self, other: Matrix) -> Matrix {
        -other.sub(self)
    }
}

impl<'a> Sub<&'a Matrix> for f64 {
    type Output = Matrix;

    fn sub(self, other: &'a Matrix) -> Self::Output {
        -other.sub(self)
    }
}

impl Sub<Matrix> for i32 {
    type Output = Matrix;

    fn sub(self, other: Matrix) -> Matrix {
        -other.sub(self)
    }
}

impl<'a> Sub<&'a Matrix> for i32 {
    type Output = Matrix;

    fn sub(self, other: &'a Matrix) -> Self::Output {
        -other.sub(self)
    }
}

impl Sub<Matrix> for usize {
    type Output = Matrix;

    fn sub(self, other: Matrix) -> Matrix {
        -other.sub(self as f64)
    }
}

impl<'a> Sub<&'a Matrix> for usize {
    type Output = Matrix;

    fn sub(self, other: &'a Matrix) -> Self::Output {
        -other.sub(self as f64)
    }
}

// =============================================================================
// Multiplication for Matrix
// =============================================================================
/// Element-wise multiplication between Matrix vs f64
impl Mul<f64> for Matrix {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &self.data;
                let mut y = vec![0f64; x.len()];
                let a_f64 = other;
                let n_i32 = x.len() as i32;

                unsafe {
                    daxpy(n_i32, a_f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => self.fmap(|x| x * other),
        }
    }
}

impl Mul<i64> for Matrix {
    type Output = Self;

    fn mul(self, other: i64) -> Self {
        self.mul(other as f64)
    }
}

impl Mul<i32> for Matrix {
    type Output = Self;

    fn mul(self, other: i32) -> Self {
        self.mul(other as f64)
    }
}

impl Mul<usize> for Matrix {
    type Output = Self;

    fn mul(self, other: usize) -> Self {
        self.mul(other as f64)
    }
}

// impl<'a> Mul<i64> for &'a Matrix {
//     type Output = Matrix;
//
//     fn mul(self, other: i64) -> Self::Output {
//         self.mul(other as f64)
//     }
// }
//
// impl<'a> Mul<i32> for &'a Matrix {
//     type Output = Matrix;
//
//     fn mul(self, other: i32) -> Self::Output {
//         self.mul(other as f64)
//     }
// }
//
// impl<'a> Mul<usize> for &'a Matrix {
//     type Output = Matrix;
//
//     fn mul(self, other: usize) -> Self::Output {
//         self.mul(other as f64)
//     }
// }

impl Mul<Matrix> for f64 {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        other.mul(self)
    }
}

impl Mul<Matrix> for i64 {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        other.mul(self as f64)
    }
}

impl Mul<Matrix> for i32 {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        other.mul(self)
    }
}

impl Mul<Matrix> for usize {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        other.mul(self as f64)
    }
}

impl<'a> Mul<&'a Matrix> for f64 {
    type Output = Matrix;

    fn mul(self, other: &'a Matrix) -> Matrix {
        other.mul_scalar(self)
    }
}

impl<'a> Mul<&'a Matrix> for i64 {
    type Output = Matrix;

    fn mul(self, other: &'a Matrix) -> Matrix {
        other.mul_scalar(self as f64)
    }
}

impl<'a> Mul<&'a Matrix> for i32 {
    type Output = Matrix;

    fn mul(self, other: &'a Matrix) -> Matrix {
        other.mul_scalar(self as f64)
    }
}

impl<'a> Mul<&'a Matrix> for usize {
    type Output = Matrix;

    fn mul(self, other: &'a Matrix) -> Matrix {
        other.mul_scalar(self as f64)
    }
}

/// Matrix Multiplication
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;4;1, 2, 2, Row);
///     let b = matrix!(1;4;1, 2, 2, Col);
///     assert_eq!(a * b, matrix(c!(5, 11, 11, 25), 2, 2, Row));
/// }
/// ```
impl Mul<Matrix> for Matrix {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        match () {
            #[cfg(feature = "O3")]
            () => blas_mul(&self, &other),
            _ => matmul(&self, &other),
        }
    }
}

impl<'a, 'b> Mul<&'b Matrix> for &'a Matrix {
    type Output = Matrix;

    fn mul(self, other: &'b Matrix) -> Self::Output {
        match () {
            #[cfg(feature = "O3")]
            () => blas_mul(self, other),
            _ => matmul(self, other),
        }
    }
}

#[allow(non_snake_case)]
impl Mul<Vec<f64>> for Matrix {
    type Output = Vec<f64>;

    fn mul(self, other: Vec<f64>) -> Self::Output {
        self.apply(&other)
    }
}

#[allow(non_snake_case)]
impl<'a, 'b> Mul<&'b Vec<f64>> for &'a Matrix {
    type Output = Vec<f64>;

    fn mul(self, other: &'b Vec<f64>) -> Self::Output {
        self.apply(other)
    }
}

/// Matrix multiplication for `Vec<f64>` vs `Matrix`
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;4;1, 2, 2, Row);
///     let v = c!(1,2);
///     assert_eq!(v * a, c!(7, 10));
/// }
/// ```
impl Mul<Matrix> for Vec<f64> {
    type Output = Vec<f64>;

    fn mul(self, other: Matrix) -> Self::Output {
        assert_eq!(self.len(), other.row);
        let mut c = vec![0f64; other.col];
        gevm(1f64, &self, &other, 0f64, &mut c);
        c
    }
}

impl<'a, 'b> Mul<&'b Matrix> for &'a Vec<f64> {
    type Output = Vec<f64>;

    fn mul(self, other: &'b Matrix) -> Self::Output {
        assert_eq!(self.len(), other.row);
        let mut c = vec![0f64; other.col];
        gevm(1f64, self, other, 0f64, &mut c);
        c
    }
}

// =============================================================================
// Standard Operation for Matrix (DIV)
// =============================================================================
/// Element-wise division between Matrix vs f64
impl Div<f64> for Matrix {
    type Output = Self;

    fn div(self, other: f64) -> Self {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &self.data;
                let mut y = vec![0f64; x.len()];
                let a_f64 = other;
                let n_i32 = x.len() as i32;

                unsafe {
                    daxpy(n_i32, 1f64 / a_f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => self.fmap(|x| x / other),
        }
    }
}

impl Div<i64> for Matrix {
    type Output = Self;

    fn div(self, other: i64) -> Self {
        self.div(other as f64)
    }
}

impl Div<i32> for Matrix {
    type Output = Self;

    fn div(self, other: i32) -> Self {
        self.div(other as f64)
    }
}

impl Div<usize> for Matrix {
    type Output = Self;

    fn div(self, other: usize) -> Self {
        self.div(other as f64)
    }
}

impl<'a> Div<f64> for &'a Matrix {
    type Output = Matrix;

    fn div(self, other: f64) -> Self::Output {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &self.data;
                let mut y = vec![0f64; x.len()];
                let a_f64 = other;
                let n_i32 = x.len() as i32;

                unsafe {
                    daxpy(n_i32, 1f64 / a_f64, x, 1, &mut y, 1);
                }
                matrix(y, self.row, self.col, self.shape)
            }
            _ => self.fmap(|x| x / other),
        }
    }
}

impl<'a> Div<i64> for &'a Matrix {
    type Output = Matrix;

    fn div(self, other: i64) -> Self::Output {
        self.div(other as f64)
    }
}

impl<'a> Div<i32> for &'a Matrix {
    type Output = Matrix;

    fn div(self, other: i32) -> Self::Output {
        self.div(other as f64)
    }
}

impl<'a> Div<usize> for &'a Matrix {
    type Output = Matrix;

    fn div(self, other: usize) -> Self::Output {
        self.div(other as f64)
    }
}

/// Index for Matrix
///
/// `(usize, usize) -> f64`
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = matrix(vec![1,2,3,4],2,2,Row);
/// assert_eq!(a[(0,1)], 2f64);
/// ```
impl Index<(usize, usize)> for Matrix {
    type Output = f64;

    fn index(&self, pair: (usize, usize)) -> &f64 {
        let p = self.ptr();
        let i = pair.0;
        let j = pair.1;
        assert!(i < self.row && j < self.col, "Index out of range");
        match self.shape {
            Row => unsafe { &*p.add(i * self.col + j) },
            Col => unsafe { &*p.add(i + j * self.row) },
        }
    }
}

/// IndexMut for Matrix (Assign)
///
/// `(usize, usize) -> f64`
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let mut a = matrix!(1;4;1, 2, 2, Row);
///     a[(1,1)] = 10.0;
///     assert_eq!(a, matrix(c!(1,2,3,10), 2, 2, Row));
/// }
/// ```
impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, pair: (usize, usize)) -> &mut f64 {
        let i = pair.0;
        let j = pair.1;
        let r = self.row;
        let c = self.col;
        assert!(i < self.row && j < self.col, "Index out of range");
        let p = self.mut_ptr();
        match self.shape {
            Row => {
                let idx = i * c + j;
                unsafe { &mut *p.add(idx) }
            }
            Col => {
                let idx = i + j * r;
                unsafe { &mut *p.add(idx) }
            }
        }
    }
}

// =============================================================================
// Functional Programming Tools (Hand-written)
// =============================================================================

impl FPMatrix for Matrix {
    fn take_row(&self, n: usize) -> Self {
        if n >= self.row {
            return self.clone();
        }
        match self.shape {
            Row => {
                let new_data = self
                    .data
                    .clone()
                    .into_iter()
                    .take(n * self.col)
                    .collect::<Vec<f64>>();
                matrix(new_data, n, self.col, Row)
            }
            Col => {
                let mut temp_data: Vec<f64> = Vec::new();
                for i in 0..n {
                    temp_data.extend(self.row(i));
                }
                matrix(temp_data, n, self.col, Row)
            }
        }
    }

    fn take_col(&self, n: usize) -> Self {
        if n >= self.col {
            return self.clone();
        }
        match self.shape {
            Col => {
                let new_data = self
                    .data
                    .clone()
                    .into_iter()
                    .take(n * self.row)
                    .collect::<Vec<f64>>();
                matrix(new_data, self.row, n, Col)
            }
            Row => {
                let mut temp_data: Vec<f64> = Vec::new();
                for i in 0..n {
                    temp_data.extend(self.col(i));
                }
                matrix(temp_data, self.row, n, Col)
            }
        }
    }

    fn skip_row(&self, n: usize) -> Self {
        assert!(n < self.row, "Skip range is larger than row of matrix");

        let mut temp_data: Vec<f64> = Vec::new();
        for i in n..self.row {
            temp_data.extend(self.row(i));
        }
        matrix(temp_data, self.row - n, self.col, Row)
    }

    fn skip_col(&self, n: usize) -> Self {
        assert!(n < self.col, "Skip range is larger than col of matrix");

        let mut temp_data: Vec<f64> = Vec::new();
        for i in n..self.col {
            temp_data.extend(self.col(i));
        }
        matrix(temp_data, self.row, self.col - n, Col)
    }

    fn fmap<F>(&self, f: F) -> Matrix
    where
        F: Fn(f64) -> f64,
    {
        let result = self.data.iter().map(|x| f(*x)).collect::<Vec<f64>>();
        matrix(result, self.row, self.col, self.shape)
    }

    /// Column map
    ///
    /// # Example
    /// ```rust
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let x = ml_matrix("1 2;3 4;5 6");
    ///     let y = x.col_map(|c| c.fmap(|t| t - c.mean()));
    ///
    ///     assert_eq!(y, ml_matrix("-2 -2;0 0;2 2"));
    /// }
    /// ```
    fn col_map<F>(&self, f: F) -> Matrix
    where
        F: Fn(Vec<f64>) -> Vec<f64>,
    {
        let mut result = matrix(vec![0f64; self.row * self.col], self.row, self.col, Col);

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
    ///
    /// fn main() {
    ///     let x = ml_matrix("1 2 3;4 5 6");
    ///     let y = x.row_map(|r| r.fmap(|t| t - r.mean()));
    ///
    ///     assert_eq!(y, ml_matrix("-1 0 1;-1 0 1"));
    /// }
    /// ```
    fn row_map<F>(&self, f: F) -> Matrix
    where
        F: Fn(Vec<f64>) -> Vec<f64>,
    {
        let mut result = matrix(vec![0f64; self.row * self.col], self.row, self.col, Row);

        for i in 0..self.row {
            result.subs_row(i, &f(self.row(i)));
        }

        result
    }

    fn col_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<f64>) -> Vec<f64>,
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
        F: Fn(Vec<f64>) -> Vec<f64>,
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

    fn reduce<F, T>(&self, init: T, f: F) -> f64
    where
        F: Fn(f64, f64) -> f64,
        T: convert::Into<f64>,
    {
        self.data.iter().fold(init.into(), |x, y| f(x, *y))
    }

    fn zip_with<F>(&self, f: F, other: &Matrix) -> Self
    where
        F: Fn(f64, f64) -> f64,
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
            .collect::<Vec<f64>>();
        matrix(result, self.row, self.col, self.shape)
    }

    fn col_reduce<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(Vec<f64>) -> f64,
    {
        let mut v = vec![0f64; self.col];
        for i in 0..self.col {
            v[i] = f(self.col(i));
        }
        v
    }

    fn row_reduce<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(Vec<f64>) -> f64,
    {
        let mut v = vec![0f64; self.row];
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
    fn back_subs(&self, b: &Vec<f64>) -> Vec<f64>;
    fn forward_subs(&self, b: &Vec<f64>) -> Vec<f64>;
    fn lu(&self) -> PQLU;
    fn waz(&self, d_form: Form) -> Option<WAZD>;
    fn qr(&self) -> QR;
    fn svd(&self) -> SVD;
    #[cfg(feature = "O3")]
    fn cholesky(&self, uplo: UPLO) -> Matrix;
    fn rref(&self) -> Matrix;
    fn det(&self) -> f64;
    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix);
    fn inv(&self) -> Matrix;
    fn pseudo_inv(&self) -> Matrix;
    fn solve(&self, b: &Vec<f64>, sk: SolveKind) -> Vec<f64>;
    fn solve_mat(&self, m: &Matrix, sk: SolveKind) -> Matrix;
    fn is_symmetric(&self) -> bool;
}

pub fn diag(n: usize) -> Matrix {
    let mut v: Vec<f64> = vec![0f64; n * n];
    for i in 0..n {
        let idx = i * (n + 1);
        v[idx] = 1f64;
    }
    matrix(v, n, n, Row)
}

/// Data structure for Complete Pivoting LU decomposition
///
/// # Usage
/// ```rust
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = ml_matrix("1 2;3 4");
/// let pqlu = a.lu();
/// let (p, q, l, u) = pqlu.extract();
/// // p, q are permutations
/// // l, u are matrices
/// l.print(); // lower triangular
/// u.print(); // upper triangular
/// ```
#[derive(Debug, Clone)]
pub struct PQLU {
    pub p: Vec<usize>,
    pub q: Vec<usize>,
    pub l: Matrix,
    pub u: Matrix,
}

impl PQLU {
    pub fn extract(&self) -> (Vec<usize>, Vec<usize>, Matrix, Matrix) {
        (
            self.p.clone(),
            self.q.clone(),
            self.l.clone(),
            self.u.clone(),
        )
    }

    pub fn det(&self) -> f64 {
        // sgn of perms
        let mut sgn_p = 1f64;
        let mut sgn_q = 1f64;
        for (i, &j) in self.p.iter().enumerate() {
            if i != j {
                sgn_p *= -1f64;
            }
        }
        for (i, &j) in self.q.iter().enumerate() {
            if i != j {
                sgn_q *= -1f64;
            }
        }

        self.u.diag().reduce(1f64, |x, y| x * y) * sgn_p * sgn_q
    }

    pub fn inv(&self) -> Matrix {
        let (p, q, l, u) = self.extract();
        let mut m = inv_u(u) * inv_l(l);
        // Q = Q1 Q2 Q3 ..
        for (idx1, idx2) in q.into_iter().enumerate().rev() {
            unsafe {
                m.swap(idx1, idx2, Row);
            }
        }
        // P = Pn-1 .. P3 P2 P1
        for (idx1, idx2) in p.into_iter().enumerate().rev() {
            unsafe {
                m.swap(idx1, idx2, Col);
            }
        }
        m
    }
}

#[derive(Debug, Clone)]
pub struct WAZD {
    pub w: Matrix,
    pub z: Matrix,
    pub d: Matrix,
}

#[derive(Debug, Copy, Clone)]
pub enum Form {
    Diagonal,
    Identity,
}

#[derive(Debug, Clone)]
pub struct QR {
    pub q: Matrix,
    pub r: Matrix,
}

impl QR {
    pub fn q(&self) -> &Matrix {
        &self.q
    }

    pub fn r(&self) -> &Matrix {
        &self.r
    }
}

#[derive(Debug, Clone)]
pub struct SVD {
    pub s: Vec<f64>,
    pub u: Matrix,
    pub vt: Matrix,
}

impl SVD {
    pub fn u(&self) -> &Matrix {
        &self.u
    }

    pub fn vt(&self) -> &Matrix {
        &self.vt
    }

    pub fn s_mat(&self) -> Matrix {
        let mut mat = zeros(self.u.col, self.vt.row);
        for i in 0 .. mat.row.min(mat.col) {
            mat[(i, i)] = self.s[i];
        }
        mat
    }

    /// Generated Truncated SVD
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let x = ml_matrix("1 2 3;4 5 6");
    /// # #[cfg(feature = "O3")] {
    ///     // Full SVD
    ///     let svd = x.svd();
    ///     svd.u().print();        // m x m matrix
    ///     svd.s_mat().print();    // m x n matrix
    ///     svd.vt().print();       // n x n matrix
    ///
    ///     // Truncated SVD
    ///     let svd2 = svd.truncated();
    ///     svd2.u().print();       // m x p matrix
    ///     svd2.s_mat().print();   // p x p matrix
    ///     svd2.vt().print();      // p x n matrix
    /// # }
    /// }
    /// ```
    pub fn truncated(&self) -> Self {
        let mut s: Vec<f64> = vec![];
        let mut u = matrix::<f64>(vec![], self.u.row, 0, Col);
        let mut vt = matrix::<f64>(vec![], 0, self.vt.col, Row);
        for (i, sig) in self.s.iter().enumerate() {
            if *sig == 0f64 {
                continue;
            } else {
                s.push(*sig);
                u.add_col_mut(&self.u.col(i));
                vt.add_row_mut(&self.vt.row(i));
            }
        }

        SVD {
            s,
            u,
            vt,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub enum SolveKind {
    LU,
    WAZ,
}

impl LinearAlgebra for Matrix {
    /// Backward Substitution for Upper Triangular
    fn back_subs(&self, b: &Vec<f64>) -> Vec<f64> {
        let n = self.col;
        let mut y = vec![0f64; n];
        y[n - 1] = b[n - 1] / self[(n - 1, n - 1)];
        for i in (0..n - 1).rev() {
            let mut s = 0f64;
            for j in i + 1..n {
                s += self[(i, j)] * y[j];
            }
            y[i] = 1f64 / self[(i, i)] * (b[i] - s);
        }
        y
    }

    /// Forward substitution for Lower Triangular
    fn forward_subs(&self, b: &Vec<f64>) -> Vec<f64> {
        let n = self.col;
        let mut y = vec![0f64; n];
        y[0] = b[0] / self[(0, 0)];
        for i in 1..n {
            let mut s = 0f64;
            for j in 0..i {
                s += self[(i, j)] * y[j];
            }
            y[i] = 1f64 / self[(i, i)] * (b[i] - s);
        }
        y
    }

    /// LU Decomposition Implements (Complete Pivot)
    ///
    /// # Description
    /// It use complete pivoting LU decomposition.
    /// You can get two permutations, and LU matrices.
    ///
    /// # Caution
    /// It returns `Option<PQLU>` - You should unwrap to obtain real value.
    /// `PQLU` has four field - `p`, `q`, `l`, `u`.
    /// `p`, `q` are permutations.
    /// `l`, `u` are matrices.
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix(vec![1,2,3,4], 2, 2, Row);
    ///     let pqlu = a.lu();
    ///     let (p,q,l,u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);
    ///     assert_eq!(p, vec![1]); // swap 0 & 1 (Row)
    ///     assert_eq!(q, vec![1]); // swap 0 & 1 (Col)
    ///     assert_eq!(l, matrix(c!(1,0,0.5,1),2,2,Row));
    ///     assert_eq!(u, matrix(c!(4,3,0,-0.5),2,2,Row));
    /// }
    /// ```
    fn lu(&self) -> PQLU {
        assert_eq!(self.col, self.row);
        let n = self.row;
        let len: usize = n * n;

        let mut l = eye(n);
        let mut u = matrix(vec![0f64; len], n, n, self.shape);

        let mut temp = self.clone();
        let (p, q) = gecp(&mut temp);
        for i in 0..n {
            for j in 0..i {
                // Inverse multiplier
                l[(i, j)] = -temp[(i, j)];
            }
            for j in i..n {
                u[(i, j)] = temp[(i, j)];
            }
        }
        // Pivoting L
        for i in 0..n - 1 {
            unsafe {
                let l_i = l.col_mut(i);
                for j in i + 1..l.col - 1 {
                    let dst = p[j];
                    std::ptr::swap(l_i[j], l_i[dst]);
                }
            }
        }
        PQLU { p, q, l, u }
    }

    fn waz(&self, d_form: Form) -> Option<WAZD> {
        match d_form {
            Form::Diagonal => {
                let n = self.row;
                let mut w = eye(n);
                let mut z = eye(n);
                let mut d = eye(n);
                let mut q = vec![0f64; n];
                let mut p = vec![0f64; n];

                for i in 0..n {
                    let m_i = self.col(i);
                    let pq = w.col(i).dot(&m_i);
                    d[(i, i)] = pq;
                    if pq == 0f64 {
                        return None;
                    }
                    for j in i + 1..n {
                        q[j] = w.col(j).dot(&m_i) / pq;
                        p[j] = z.col(j).dot(&self.row(i)) / pq;
                    }
                    for j in i + 1..n {
                        for k in 0..i + 1 {
                            w[(k, j)] -= q[j] * w[(k, i)];
                            z[(k, j)] -= p[j] * z[(k, i)];
                        }
                    }
                }
                Some(WAZD { w, z, d })
            }
            Form::Identity => {
                let n = self.row;
                let mut w = eye(n);
                let mut z = eye(n);
                let mut p = zeros(n, n);
                let mut q = zeros(n, n);

                for i in 0..n {
                    let m_i = self.col(i);
                    let p_ii = w.col(i).dot(&m_i);
                    p[(i, i)] = p_ii;
                    if p_ii == 0f64 {
                        return None;
                    }
                    for j in i + 1..n {
                        q[(i, j)] = w.col(j).dot(&m_i) / p_ii;
                        p[(i, j)] = z.col(j).dot(&self.row(i)) / p_ii;
                        for k in 0..j {
                            w[(k, j)] -= q[(i, j)] * w[(k, i)];
                            z[(k, j)] -= p[(i, j)] * z[(k, i)];
                        }
                    }
                    unsafe {
                        let col_ptr = z.col_mut(i);
                        col_ptr.into_iter().for_each(|x| *x /= p_ii);
                    }
                }
                Some(WAZD { w, z, d: eye(n) })
            }
        }
    }

    /// QR Decomposition
    ///
    /// Translation of [RosettaCode#Python](https://rosettacode.org/wiki/QR_decomposition#Python)
    ///
    /// # Example
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = ml_matrix("12 -51 4;6 167 -68; -4 24 -41");
    ///     let qr = a.qr();
    ///     let r = ml_matrix("-14 -21 14; 0 -175 70; 0 0 -35");
    ///     #[cfg(feature="O3")]
    ///     {
    ///         assert_eq!(r, qr.r);
    ///     }
    ///     qr.r.print();
    /// }
    /// ```
    #[allow(non_snake_case)]
    fn qr(&self) -> QR {
        match () {
            #[cfg(feature="O3")]
            () => {
                let opt_dgeqrf = lapack_dgeqrf(self);
                match opt_dgeqrf {
                    None => panic!("Serious problem in QR decomposition"),
                    Some(dgeqrf) => {
                        let q = dgeqrf.get_Q();
                        let r = dgeqrf.get_R();
                        QR {
                            q,
                            r
                        }
                    }
                }
            }
            _ => {
                let m = self.row;
                let n = self.col;

                let mut r = self.clone();
                let mut q = eye(m);
                let sub = if m == n { 1 } else { 0 };
                for i in 0..n - sub {
                    let mut H = eye(m);
                    let hh = gen_householder(&self.col(i).skip(i));
                    for j in i..m {
                        for k in i..m {
                            H[(j, k)] = hh[(j - i, k - i)];
                        }
                    }
                    q = &q * &H;
                    r = &H * &r;
                }

                QR { q, r }
            }
        }
    }

    /// Singular Value Decomposition
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = ml_matrix("3 2 2;2 3 -2");
    ///     #[cfg(feature="O3")]
    ///     {
    ///         let svd = a.svd();
    ///         assert!(eq_vec(&vec![5f64, 3f64], &svd.s, 1e-7));
    ///     }
    ///     a.print();
    /// }
    /// ```
    fn svd(&self) -> SVD {
        match () {
            #[cfg(feature="O3")]
            () => {
                let opt_dgesvd = lapack_dgesvd(self);
                match opt_dgesvd {
                    None => panic!("Illegal value in LAPACK SVD"),
                    Some(dgesvd) => match dgesvd.status {
                        SVD_STATUS::Diverge(i) => {
                            panic!("Divergent occurs in SVD - {} iterations", i)
                        }
                        SVD_STATUS::Success => {
                            SVD {
                                s: dgesvd.s,
                                u: dgesvd.u,
                                vt: dgesvd.vt,
                            }
                        }
                    }
                }
            }
            _ => {
                unimplemented!()
            }
        }
    }

    /// Cholesky Decomposition
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = ml_matrix("1 2;2 5");
    ///     #[cfg(feature = "O3")]
    ///     {
    ///         let u = a.cholesky(Upper);
    ///         let l = a.cholesky(Lower);
    ///
    ///         assert_eq!(u, ml_matrix("1 2;0 1"));
    ///         assert_eq!(l, ml_matrix("1 0;2 1"));
    ///     }
    ///     a.print();
    /// }
    /// ```
    #[cfg(feature = "O3")]
    fn cholesky(&self, uplo: UPLO) -> Matrix {
        match () {
            #[cfg(feature = "O3")]
            () => {
                if !self.is_symmetric() {
                    panic!("Cholesky Error: Matrix is not symmetric!");
                }
                let dpotrf = lapack_dpotrf(self, uplo);
                match dpotrf {
                    None => panic!("Cholesky Error: Not symmetric or not positive definite."),
                    Some(x) => {
                        match x.status {
                            POSITIVE_STATUS::Failed(i) => panic!("Cholesky Error: the leading minor of order {} is not positive definite", i),
                            POSITIVE_STATUS::Success => {
                                match uplo {
                                    UPLO::Upper => x.get_U().unwrap(),
                                    UPLO::Lower => x.get_L().unwrap()
                                }
                            }
                        }
                    }
                }
            }
            _ => {
                unimplemented!()
            }
        }
    }

    /// Reduced Row Echelon Form
    ///
    /// Implementation of [RosettaCode](https://rosettacode.org/wiki/Reduced_row_echelon_form)
    fn rref(&self) -> Matrix {
        let mut lead = 0usize;
        let mut result = self.clone();
        'outer: for r in 0..self.row {
            if self.col <= lead {
                break;
            }
            let mut i = r;
            while result[(i, lead)] == 0f64 {
                i += 1;
                if self.row == i {
                    i = r;
                    lead += 1;
                    if self.col == lead {
                        break 'outer;
                    }
                }
            }
            unsafe {
                result.swap(i, r, Row);
            }
            let tmp = result[(r, lead)];
            if tmp != 0f64 {
                unsafe {
                    result.row_mut(r).iter_mut().for_each(|t| *(*t) /= tmp);
                }
            }
            for j in 0..result.row {
                if j != r {
                    let tmp1 = result.row(r).mul_scalar(result[(j, lead)]);
                    let tmp2 = result.row(j).sub_vec(&tmp1);
                    result.subs_row(j, &tmp2);
                }
            }
            lead += 1;
        }
        result
    }

    /// Determinant
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix!(1;4;1, 2, 2, Row);
    ///     assert_eq!(a.det(), -2f64);
    /// }
    /// ```
    fn det(&self) -> f64 {
        assert_eq!(self.row, self.col);
        match () {
            #[cfg(feature = "O3")]
            () => {
                let opt_dgrf = lapack_dgetrf(self);
                match opt_dgrf {
                    None => NAN,
                    Some(dgrf) => match dgrf.status {
                        LAPACK_STATUS::Singular => 0f64,
                        LAPACK_STATUS::NonSingular => {
                            let mat = &dgrf.fact_mat;
                            let ipiv = &dgrf.ipiv;
                            let n = mat.col;
                            let mut sgn = 1i32;
                            let mut d = 1f64;
                            for i in 0..n {
                                d *= mat[(i, i)];
                            }
                            for i in 0..ipiv.len() {
                                if ipiv[i] - 1 != i as i32 {
                                    sgn *= -1;
                                }
                            }
                            (sgn as f64) * d
                        }
                    },
                }
            }
            _ => self.lu().det(),
        }
    }

    /// Block Partition
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix!(1;16;1, 4, 4, Row);
    ///     let (m1, m2, m3, m4) = a.block();
    ///     assert_eq!(m1, matrix(c!(1,2,5,6), 2, 2, Row));
    ///     assert_eq!(m2, matrix(c!(3,4,7,8), 2, 2, Row));
    ///     assert_eq!(m3, matrix(c!(9,10,13,14), 2, 2, Row));
    ///     assert_eq!(m4, matrix(c!(11,12,15,16), 2, 2, Row));
    ///
    ///     let b = matrix!(1;16;1, 4, 4, Col);
    ///     let (m1, m2, m3, m4) = b.block();
    ///     assert_eq!(m1, matrix(c!(1,2,5,6), 2, 2, Col));
    ///     assert_eq!(m3, matrix(c!(3,4,7,8), 2, 2, Col));
    ///     assert_eq!(m2, matrix(c!(9,10,13,14), 2, 2, Col));
    ///     assert_eq!(m4, matrix(c!(11,12,15,16), 2, 2, Col));
    /// }
    /// ```
    fn block(&self) -> (Self, Self, Self, Self) {
        let r = self.row;
        let c = self.col;
        let l_r = self.row / 2;
        let l_c = self.col / 2;
        let r_l = r - l_r;
        let c_l = c - l_c;

        let mut m1 = matrix(vec![0f64; l_r * l_c], l_r, l_c, self.shape);
        let mut m2 = matrix(vec![0f64; l_r * c_l], l_r, c_l, self.shape);
        let mut m3 = matrix(vec![0f64; r_l * l_c], r_l, l_c, self.shape);
        let mut m4 = matrix(vec![0f64; r_l * c_l], r_l, c_l, self.shape);

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

    /// Inverse of Matrix
    ///
    /// # Caution
    ///
    /// `inv` function returns `Option<Matrix>`
    /// Thus, you should use pattern matching or `unwrap` to obtain inverse.
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     // Non-singular
    ///     let a = matrix!(1;4;1, 2, 2, Row);
    ///     assert_eq!(a.inv(), matrix(c!(-2,1,1.5,-0.5),2,2,Row));
    ///
    ///     // Singular
    ///     //let b = matrix!(1;9;1, 3, 3, Row);
    ///     //let c = b.inv(); // Can compile but..
    ///     //let d = b.det();
    ///     //assert_eq!(d, 0f64);
    /// }
    /// ```
    fn inv(&self) -> Self {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let opt_dgrf = lapack_dgetrf(self);
                match opt_dgrf {
                    None => panic!("Singular matrix (opt_dgrf)"),
                    Some(dgrf) => match dgrf.status {
                        LAPACK_STATUS::NonSingular => lapack_dgetri(&dgrf).unwrap(),
                        LAPACK_STATUS::Singular => panic!("Singular matrix (LAPACK_STATUS Singular)"),
                    },
                }
            }
            _ => self.lu().inv(),
            // _ => {
            //     match self.lu() {
            //         None => None,
            //         Some(lu) => Some(lu.inv())
            //     }
            // }
        }
    }

    /// Moore-Penrose Pseudo inverse
    ///
    /// # Description
    /// `$X^\dagger = (X^T X)^{-1} X^T$`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = matrix!(1;4;1, 2, 2, Row);
    ///     let inv_a = a.inv();
    ///     let pse_a = a.pseudo_inv();
    ///
    ///     assert_eq!(inv_a, pse_a); // Nearly equal
    /// }
    /// ```
    fn pseudo_inv(&self) -> Self {
        match () {
            #[cfg(feature="O3")]
            () => {
                let svd = self.svd();
                let row = svd.vt.row;
                let col = svd.u.row;
                let mut sp = zeros(row, col);
                for i in 0 .. row.min(col) {
                    sp[(i, i)] = 1f64 / svd.s[i];
                }
                svd.vt.t() * sp * svd.u.t()
            }
            _ => {
                let xt = self.t();
                let xtx = &xt * self;
                xtx.inv() * xt
            }
        }
    }

    /// Solve with Vector
    ///
    /// # Solve options
    ///
    /// * LU: Gaussian elimination with Complete pivoting LU (GECP)
    /// * WAZ: Solve with WAZ decomposition
    ///
    /// # Reference
    ///
    /// * Biswa Nath Datta, *Numerical Linear Algebra and Applications, Second Edition*
    /// * Ke Chen, *Matrix Preconditioning Techniques and Applications*, Cambridge Monographs on Applied and Computational Mathematics
    fn solve(&self, b: &Vec<f64>, sk: SolveKind) -> Vec<f64> {
        match sk {
            #[cfg(feature = "O3")]
            SolveKind::LU => {
                let opt_dgrf = lapack_dgetrf(self);
                match opt_dgrf {
                    None => panic!("Try solve for Singluar matrix"),
                    Some(dgrf) => match dgrf.status {
                        LAPACK_STATUS::Singular => panic!("Try solve for Singluar matrix"),
                        LAPACK_STATUS::NonSingular => {
                            lapack_dgetrs(&dgrf, &(b.into())).unwrap().into()
                        }
                    },
                }
            }
            #[cfg(not(feature = "O3"))]
            SolveKind::LU => {
                let lu = self.lu();
                let (p, q, l, u) = lu.extract();
                let mut v = b.clone();
                v.swap_with_perm(&p.into_iter().enumerate().collect());
                let z = l.forward_subs(&v);
                let mut y = u.back_subs(&z);
                y.swap_with_perm(&q.into_iter().enumerate().rev().collect());
                y
            }
            SolveKind::WAZ => {
                let wazd = match self.waz(Form::Identity) {
                    None => panic!("Can't solve by WAZ with Singular matrix!"),
                    Some(obj) => obj,
                };
                let x = &wazd.w.t() * b;
                let x = &wazd.z * &x;
                x
            }
        }
    }

    fn solve_mat(&self, m: &Matrix, sk: SolveKind) -> Matrix {
        match sk {
            #[cfg(feature = "O3")]
            SolveKind::LU => {
                let opt_dgrf = lapack_dgetrf(self);
                match opt_dgrf {
                    None => panic!("Try solve for Singluar matrix"),
                    Some(dgrf) => match dgrf.status {
                        LAPACK_STATUS::Singular => panic!("Try solve for Singluar matrix"),
                        LAPACK_STATUS::NonSingular => lapack_dgetrs(&dgrf, m).unwrap(),
                    },
                }
            }
            #[cfg(not(feature = "O3"))]
            SolveKind::LU => {
                let lu = self.lu();
                let (p, q, l, u) = lu.extract();
                let mut x = matrix(vec![0f64; self.col * m.col], self.col, m.col, Col);
                for i in 0..m.col {
                    let mut v = m.col(i).clone();
                    for (r, &s) in p.iter().enumerate() {
                        v.swap(r, s);
                    }
                    let z = l.forward_subs(&v);
                    let mut y = u.back_subs(&z);
                    for (r, &s) in q.iter().enumerate() {
                        y.swap(r, s);
                    }
                    unsafe {
                        let mut c = x.col_mut(i);
                        copy_vec_ptr(&mut c, &y);
                    }
                }
                x
            }
            SolveKind::WAZ => {
                let wazd = match self.waz(Form::Identity) {
                    None => panic!("Try solve for Singular matrix"),
                    Some(obj) => obj,
                };
                let x = &wazd.w.t() * m;
                let x = &wazd.z * &x;
                x
            }
        }
    }

    fn is_symmetric(&self) -> bool {
        if self.row != self.col {
            return false;
        }

        for i in 0 .. self.row {
            for j in i .. self.col {
                if !nearly_eq(self[(i,j)], self[(j,i)]) {
                    return false;
                }
            }
        }
        true
    }
}

#[allow(non_snake_case)]
pub fn solve(A: &Matrix, b: &Matrix, sk: SolveKind) -> Matrix {
    A.solve_mat(b, sk)
}

impl MutMatrix for Matrix {
    unsafe fn col_mut(&mut self, idx: usize) -> Vec<*mut f64> {
        assert!(idx < self.col, "Index out of range");
        match self.shape {
            Shape::Col => {
                let mut v: Vec<*mut f64> = Vec::with_capacity(self.row);
                v.set_len(self.row);
                let start_idx = idx * self.row;
                let p = self.mut_ptr();
                for (i, j) in (start_idx..start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Row => {
                let mut v: Vec<*mut f64> = Vec::with_capacity(self.row);
                v.set_len(self.row);
                let p = self.mut_ptr();
                for i in 0..v.len() {
                    v[i] = p.add(idx + i * self.col);
                }
                v
            }
        }
    }

    unsafe fn row_mut(&mut self, idx: usize) -> Vec<*mut f64> {
        assert!(idx < self.row, "Index out of range");
        match self.shape {
            Shape::Row => {
                let mut v: Vec<*mut f64> = Vec::with_capacity(self.col);
                v.set_len(self.col);
                let start_idx = idx * self.col;
                let p = self.mut_ptr();
                for (i, j) in (start_idx..start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Col => {
                let mut v: Vec<*mut f64> = Vec::with_capacity(self.col);
                v.set_len(self.col);
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
            Shape::Col => swap_vec_ptr(&mut self.col_mut(idx1), &mut self.col_mut(idx2)),
            Shape::Row => swap_vec_ptr(&mut self.row_mut(idx1), &mut self.row_mut(idx2)),
        }
    }

    unsafe fn swap_with_perm(&mut self, p: &Vec<(usize, usize)>, shape: Shape) {
        for (i, j) in p.iter() {
            self.swap(*i, *j, shape)
        }
    }
}

// =============================================================================
// Back-end Utils
// =============================================================================

/// Combine separated matrix to one matrix
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;16;1, 4, 4, Row);
///     let (m1, m2, m3, m4) = a.block();
///     let m = combine(m1,m2,m3,m4);
///     assert_eq!(m, a);
///
///     let b = matrix!(1;16;1, 4, 4, Col);
///     let (n1, n2, n3, n4) = b.block();
///     let n = combine(n1,n2,n3,n4);
///     assert_eq!(n, b);
/// }
/// ```
pub fn combine(m1: Matrix, m2: Matrix, m3: Matrix, m4: Matrix) -> Matrix {
    let l_r = m1.row;
    let l_c = m1.col;
    let c_l = m2.col;
    let r_l = m3.row;

    let r = l_r + r_l;
    let c = l_c + c_l;

    let mut m = matrix(vec![0f64; r * c], r, c, m1.shape);

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
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix(c!(1,0,2,1), 2, 2, Row);
///     assert_eq!(inv_l(a), matrix(c!(1,0,-2,1), 2, 2, Row));
///
///     let b = matrix(c!(1,0,0,2,1,0,4,3,1), 3, 3, Row);
///     assert_eq!(inv_l(b), matrix(c!(1,0,0,-2,1,0,2,-3,1), 3, 3, Row));
/// }
/// ```
pub fn inv_l(l: Matrix) -> Matrix {
    let mut m = l.clone();

    match l.row {
        1 => l,
        2 => {
            m[(1, 0)] = -m[(1, 0)];
            m
        }
        _ => {
            let (l1, l2, l3, l4) = l.block();

            let m1 = inv_l(l1);
            let m2 = l2;
            let m4 = inv_l(l4);
            let m3 = -(&(&m4 * &l3) * &m1);

            combine(m1, m2, m3, m4)
        }
    }
}

/// Inverse of upper triangular matrix
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let u = matrix(c!(2,2,0,1), 2, 2, Row);
///     assert_eq!(inv_u(u), matrix(c!(0.5,-1,0,1), 2, 2, Row));
/// }
/// ```
pub fn inv_u(u: Matrix) -> Matrix {
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
            let m1 = inv_u(u1);
            let m3 = u3;
            let m4 = inv_u(u4);
            let m2 = -(m1.clone() * u2 * m4.clone());

            combine(m1, m2, m3, m4)
        }
    }
}

/// Matrix multiply back-ends
fn matmul(a: &Matrix, b: &Matrix) -> Matrix {
    assert_eq!(a.col, b.row);
    let mut c = matrix(vec![0f64; a.row * b.col], a.row, b.col, a.shape);
    gemm(1f64, a, b, 0f64, &mut c);
    c
}

/// GEMM wrapper for Matrixmultiply
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::prelude::*;
///
/// fn main() {
///     let a = ml_matrix("1 2 3;4 5 6");
///     let b = ml_matrix("1 2;3 4;5 6");
///     let mut c1 = zeros(2, 2);
///     let mut c2 = matrix(vec![1f64; 9], 3, 3, Col);
///
///     gemm(1f64, &a, &b, 0f64, &mut c1);
///     gemm(1f64, &b, &a, 2f64, &mut c2);
///
///     assert_eq!(c1, ml_matrix("22 28; 49 64"));
///     assert_eq!(c2, ml_matrix("11 14 17;21 28 35;31 42 53"));
/// }
/// ```
pub fn gemm(alpha: f64, a: &Matrix, b: &Matrix, beta: f64, c: &mut Matrix) {
    let m = a.row;
    let k = a.col;
    let n = b.col;
    let (rsa, csa) = match a.shape {
        Row => (a.col as isize, 1isize),
        Col => (1isize, a.row as isize),
    };
    let (rsb, csb) = match b.shape {
        Row => (b.col as isize, 1isize),
        Col => (1isize, b.row as isize),
    };
    let (rsc, csc) = match c.shape {
        Row => (c.col as isize, 1isize),
        Col => (1isize, c.row as isize),
    };

    unsafe {
        matrixmultiply::dgemm(
            m,
            k,
            n,
            alpha,
            a.ptr(),
            rsa,
            csa,
            b.ptr(),
            rsb,
            csb,
            beta,
            c.mut_ptr(),
            rsc,
            csc,
        )
    }
}

/// General Matrix-Vector multiplication
///
/// # Example
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = ml_matrix("1 2 3; 4 5 6");
///     let b = c!(1, 2, 3);
///     let mut c = vec![0f64; 2];
///     gemv(1f64, &a, &b, 0f64, &mut c);
///     assert_eq!(c, c!(14, 32));
/// }
/// ```
pub fn gemv(alpha: f64, a: &Matrix, b: &Vec<f64>, beta: f64, c: &mut Vec<f64>) {
    let m = a.row;
    let k = a.col;
    let n = 1usize;
    let (rsa, csa) = match a.shape {
        Row => (a.col as isize, 1isize),
        Col => (1isize, a.row as isize),
    };
    let (rsb, csb) = (1isize, 1isize);
    let (rsc, csc) = (1isize, 1isize);

    unsafe {
        matrixmultiply::dgemm(
            m,
            k,
            n,
            alpha,
            a.ptr(),
            rsa,
            csa,
            b.as_ptr(),
            rsb,
            csb,
            beta,
            c.as_mut_ptr(),
            rsc,
            csc,
        )
    }
}

/// General Vector-Matrix multiplication
///
/// # Example
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = c!(1, 2);
///     let b = ml_matrix("1 2 3; 4 5 6");
///     let mut c = vec![0f64; 3];
///     gevm(1f64, &a, &b, 0f64, &mut c);
///     assert_eq!(c, c!(9, 12, 15));
/// }
/// ```
pub fn gevm(alpha: f64, a: &Vec<f64>, b: &Matrix, beta: f64, c: &mut Vec<f64>) {
    let m = 1usize;
    let k = a.len();
    let n = b.col;
    let (rsa, csa) = (1isize, 1isize);
    let (rsb, csb) = match b.shape {
        Row => (b.col as isize, 1isize),
        Col => (1isize, b.row as isize),
    };
    let (rsc, csc) = (1isize, 1isize);

    unsafe {
        matrixmultiply::dgemm(
            m,
            k,
            n,
            alpha,
            a.as_ptr(),
            rsa,
            csa,
            b.ptr(),
            rsb,
            csb,
            beta,
            c.as_mut_ptr(),
            rsc,
            csc,
        )
    }
}

//fn matmul(a: &Matrix, b: &Matrix) -> Matrix {
//    match (a.row, a.col) {
//        (p, q) if p <= 100 && q <= 100 => {
//            let r_self = a.row;
//            let c_self = a.col;
//            let new_other = b;
//            let r_other = new_other.row;
//            let c_other = new_other.col;
//
//            assert_eq!(c_self, r_other);
//
//            let r_new = r_self;
//            let c_new = c_other;
//
//            let mut result = matrix(vec![0f64; r_new * c_new], r_new, c_new, a.shape);
//
//            for i in 0..r_new {
//                for j in 0..c_new {
//                    let mut s = 0f64;
//                    for k in 0..c_self {
//                        s += a[(i, k)] * new_other[(k, j)];
//                    }
//                    result[(i, j)] = s;
//                }
//            }
//            result
//        }
//        _ => {
//            let (a1, a2, a3, a4) = a.block();
//            let (b1, b2, b3, b4) = b.block();
//
//            let m1 = matmul(&a1, &b1) + matmul(&a2, &b3);
//            let m2 = matmul(&a1, &b2) + matmul(&a2, &b4);
//            let m3 = matmul(&a3, &b1) + matmul(&a4, &b3);
//            let m4 = matmul(&a3, &b2) + matmul(&a4, &b4);
//
//            combine(m1, m2, m3, m4)
//        }
//    }
//}

// =============================================================================
// BLAS & LAPACK Area
// =============================================================================

/// Matrix multiplication with BLAS
///
/// * m1: m x k matrix
/// * m2: k x n matrix
/// * result: m x n matrix
#[cfg(feature = "O3")]
pub fn blas_mul(m1: &Matrix, m2: &Matrix) -> Matrix {
    let m = m1.row;
    let k = m1.col;
    assert_eq!(k, m2.row);
    let n = m2.col;

    let m_i32 = m as i32;
    let n_i32 = n as i32;
    let k_i32 = k as i32;

    let a = &m1.data;
    let b = &m2.data;
    let mut c = vec![0f64; m * n];

    match (m1.shape, m2.shape) {
        (Row, Row) => {
            unsafe {
                dgemm(
                    b'N', b'N', n_i32, m_i32, k_i32, 1f64, b, n_i32, a, k_i32, 0f64, &mut c, n_i32,
                );
            }
            matrix(c, m, n, Row)
        }
        (Row, Col) => {
            unsafe {
                dgemm(
                    b'T', b'N', n_i32, m_i32, k_i32, 1f64, b, k_i32, a, k_i32, 0f64, &mut c, n_i32,
                );
            }
            matrix(c, m, n, Row)
        }
        (Col, Col) => {
            unsafe {
                // (nxk) x (kxm) = n x m
                dgemm(
                    b'N', b'N', m_i32, n_i32, k_i32, 1f64, a, m_i32, b, k_i32, 0f64, &mut c, m_i32,
                );
            }
            matrix(c, m, n, Col)
        }
        (Col, Row) => {
            unsafe {
                dgemm(
                    b'N', b'T', m_i32, n_i32, k_i32, 1f64, a, m_i32, b, n_i32, 0f64, &mut c, m_i32,
                );
            }
            matrix(c, m, n, Col)
        }
    }
}

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum LAPACK_STATUS {
    Singular,
    NonSingular,
}

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum SVD_STATUS {
    Success,
    Diverge(i32),
}

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum POSITIVE_STATUS {
    Success,
    Failed(i32),
}

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum UPLO {
    Upper,
    Lower
}

/// Temporary data structure from `dgetrf`
#[derive(Debug, Clone)]
pub struct DGETRF {
    pub fact_mat: Matrix,
    pub ipiv: Vec<i32>,
    pub status: LAPACK_STATUS,
}

/// Temporary data structure from `dgeqrf`
#[derive(Debug, Clone)]
pub struct DGEQRF {
    pub fact_mat: Matrix,
    pub tau: Vec<f64>,
    pub status: LAPACK_STATUS,
}

#[derive(Debug, Clone)]
pub struct DGESVD {
    pub s: Vec<f64>,
    pub u: Matrix,
    pub vt: Matrix,
    pub status: SVD_STATUS,
}

#[derive(Debug, Clone)]
pub struct DPOTRF {
    pub fact_mat: Matrix,
    pub uplo: UPLO,
    pub status: POSITIVE_STATUS
}

///// Temporary data structure from `dgeev`
//#[derive(Debug, Clone)]
//pub struct DGEEV {
//    pub fact_mat: Matrix,
//    pub tau: Vec<f64>,
//    pub status: LAPACK_STATUS,
//}

/// Peroxide version of `dgetrf`
#[cfg(feature = "O3")]
pub fn lapack_dgetrf(mat: &Matrix) -> Option<DGETRF> {
    let m = mat.row;
    let n = mat.col;
    let m_i32 = m as i32;
    let n_i32 = n as i32;
    // Should column major
    let mut a = match mat.shape {
        Row => mat.change_shape().data.clone(),
        Col => mat.data.clone(),
    };

    let mut info = 0i32;
    let mut ipiv: Vec<i32> = vec![0i32; std::cmp::min(m, n)];

    unsafe {
        dgetrf(m_i32, n_i32, &mut a, m_i32, &mut ipiv, &mut info);
    }

    if info < 0 {
        None
    } else if info == 0 {
        Some(DGETRF {
            fact_mat: matrix(a, m, n, Col),
            ipiv,
            status: LAPACK_STATUS::NonSingular,
        })
    } else {
        Some(DGETRF {
            fact_mat: matrix(a, m, n, Col),
            ipiv,
            status: LAPACK_STATUS::Singular,
        })
    }
}

/// Peroxide version of `dgetri`
#[cfg(feature = "O3")]
pub fn lapack_dgetri(dgrf: &DGETRF) -> Option<Matrix> {
    let mut result = dgrf.fact_mat.clone();
    let ipiv = &dgrf.ipiv;
    let (n, lda) = (result.col as i32, result.row as i32);
    let mut info = 0i32;
    let mut work = vec![0f64; result.col];

    // Workspace Query
    unsafe {
        dgetri(n, &mut result.data, lda, ipiv, &mut work, -1, &mut info);
    }

    let optimal_lwork = work[0] as usize;
    let mut optimal_work = vec![0f64; optimal_lwork];

    // Real dgetri
    unsafe {
        dgetri(
            n,
            &mut result.data,
            lda,
            ipiv,
            &mut optimal_work,
            optimal_lwork as i32,
            &mut info,
        );
    }

    if info == 0 {
        Some(result)
    } else {
        None
    }
}

/// Peroxide version of `dgetrs`
#[allow(non_snake_case)]
#[cfg(feature = "O3")]
pub fn lapack_dgetrs(dgrf: &DGETRF, b: &Matrix) -> Option<Matrix> {
    match b.shape {
        Row => lapack_dgetrs(dgrf, &b.change_shape()),
        Col => {
            let A = &dgrf.fact_mat;
            let mut b_vec = b.data.clone();
            let ipiv = &dgrf.ipiv;
            let n = A.col as i32;
            let nrhs = b.col as i32;
            let lda = A.row as i32;
            let ldb = b.row as i32;
            let mut info = 0i32;

            unsafe {
                dgetrs(
                    b'N', n, nrhs, &A.data, lda, ipiv, &mut b_vec, ldb, &mut info,
                );
            }

            if info == 0 {
                Some(matrix(b_vec, A.col, b.col, Col))
            } else {
                None
            }
        }
    }
}

/// Peroxide version of `dgeqrf`
#[allow(non_snake_case)]
#[cfg(feature = "O3")]
pub fn lapack_dgeqrf(mat: &Matrix) -> Option<DGEQRF> {
    match mat.shape {
        Row => lapack_dgeqrf(&mat.change_shape()),
        Col => {
            let m = mat.row as i32;
            let n = mat.col as i32;
            let mut A = mat.clone();
            let mut tau = vec![0f64; min(mat.row, mat.col)];
            let mut work = vec![0f64; mat.col];
            let mut info = 0i32;

            // Workspace query
            unsafe {
                dgeqrf(m, n, &mut A.data, m, &mut tau, &mut work, -1, &mut info);
            }

            let optimal_lwork = work[0] as usize;
            let mut optimal_work = vec![0f64; optimal_lwork];

            // Real dgeqrf
            unsafe {
                dgeqrf(
                    m,
                    n,
                    &mut A.data,
                    m,
                    &mut tau,
                    &mut optimal_work,
                    optimal_lwork as i32,
                    &mut info,
                );
            }

            if info == 0 {
                Some(DGEQRF {
                    fact_mat: A,
                    tau,
                    status: LAPACK_STATUS::NonSingular,
                })
            } else if info > 0 {
                Some(DGEQRF {
                    fact_mat: A,
                    tau,
                    status: LAPACK_STATUS::Singular,
                })
            } else {
                None
            }
        }
    }
}

#[allow(non_snake_case)]
#[cfg(feature = "O3")]
pub fn lapack_dgesvd(mat: &Matrix) -> Option<DGESVD> {
    match mat.shape {
        Row => lapack_dgesvd(&mat.change_shape()),
        Col => {
            let jobu = b'A';
            let jobvt = b'A';
            let m = mat.row as i32;
            let n = mat.col as i32;
            let mut A = mat.clone();
            let lda = m;
            let mut s = vec![0f64; m.min(n) as usize];
            let ldu = m;
            let mut u = vec![0f64; (ldu * m) as usize];
            let ldvt = n;
            let mut vt = vec![0f64; (ldvt * n) as usize];
            let mut work = vec![0f64; mat.col];
            let lwork = -1i32;
            let mut info = 0i32;

            // Workspace query
            unsafe {
                dgesvd(jobu, jobvt, m, n, &mut A.data, lda, &mut s, &mut u, ldu, &mut vt, ldvt, &mut work, lwork, &mut info);
            }

            let optimal_lwork = work[0] as usize;
            let mut optimal_work = vec![0f64; optimal_lwork];

            // Real dgesvd
            unsafe {
                dgesvd(
                    jobu,
                    jobvt,
                    m,
                    n,
                    &mut A.data,
                    lda,
                    &mut s,
                    &mut u,
                    ldu,
                    &mut vt,
                    ldvt,
                    &mut optimal_work,
                    optimal_lwork as i32,
                    &mut info,
                )
            }

            if info == 0 {
                Some(DGESVD {
                    s: s,
                    u: matrix(u, m as usize, m as usize, Col),
                    vt: matrix(vt, n as usize, n as usize, Col),
                    status: SVD_STATUS::Success,
                })
            } else if info < 0 {
                None
            } else {
                Some(DGESVD {
                    s: s,
                    u: matrix(u, m as usize, m as usize, Col),
                    vt: matrix(vt, n as usize, n as usize, Col),
                    status: SVD_STATUS::Diverge(info),
                })
            }
        }
    }
}

#[allow(non_snake_case)]
#[cfg(feature = "O3")]
pub fn lapack_dpotrf(mat: &Matrix, UPLO: UPLO) -> Option<DPOTRF> {
    match mat.shape {
        Row => lapack_dpotrf(&mat.change_shape(), UPLO),
        Col => {
            let lda = mat.row as i32;
            let N = mat.col as i32;
            let mut A = mat.clone();
            let mut info = 0i32;
            let uplo = match UPLO {
                UPLO::Upper => b'U',
                UPLO::Lower => b'L'
            };

            unsafe {
                dpotrf(
                    uplo,
                    N,
                    &mut A.data,
                    lda,
                    &mut info,
                )
            }

            if info == 0 {
                Some(
                    DPOTRF {
                        fact_mat: matrix(A.data, mat.row, mat.col, Col),
                        uplo: UPLO,
                        status: POSITIVE_STATUS::Success
                    }
                )
            } else if info > 0 {
                Some(
                    DPOTRF {
                        fact_mat: matrix(A.data, mat.row, mat.col, Col),
                        uplo: UPLO,
                        status: POSITIVE_STATUS::Failed(info)
                    }
                )
            } else {
                None
            }
        }
    }
}

#[allow(non_snake_case)]
#[cfg(feature = "O3")]
impl DGETRF {
    pub fn get_P(&self) -> Vec<i32> {
        self.ipiv.clone()
    }

    pub fn get_L(&self) -> Matrix {
        let mut l = self.fact_mat.clone();
        for i in 0..l.row {
            l[(i, i)] = 1f64;
            for j in i + 1..l.col {
                l[(i, j)] = 0f64;
            }
        }
        l
    }

    pub fn get_U(&self) -> Matrix {
        let mut u = self.fact_mat.clone();
        for i in 0..u.row {
            for j in 0..i {
                u[(i, j)] = 0f64;
            }
        }
        u
    }

    pub fn get_cond(&self) -> Option<f64> {
        let A = &self.fact_mat;
        let lda = A.row as i32;
        let n = A.col as i32;
        let anorm = A.norm(Norm::L1);
        let mut work = vec![0f64; 4 * A.col];
        let mut iwork = vec![0i32; A.col];
        let mut info = 0i32;
        let mut rcond = 0f64;

        unsafe {
            dgecon(
                b'1', n, &A.data, lda, anorm, &mut rcond, &mut work, &mut iwork, &mut info,
            );
        }

        if info == 0 {
            Some(rcond)
        } else {
            None
        }
    }
}

#[allow(non_snake_case)]
#[cfg(feature = "O3")]
impl DGEQRF {
    pub fn get_Q(&self) -> Matrix {
        let mut A = self.fact_mat.clone();
        let m = A.row as i32;
        let n = A.col as i32;
        let k = min(m, n);
        let lda = m;
        let tau = &self.tau;
        let mut lwork = -1i32;
        let mut work = vec![0f64; 1];
        let mut info = 0i32;

        // Optimize
        unsafe {
            dorgqr(m, n, k, &mut A.data, lda, tau, &mut work, lwork, &mut info);
        }

        lwork = work[0] as i32;
        work = vec![0f64; lwork as usize];

        // Real dorgqr
        unsafe {
            dorgqr(m, n, k, &mut A.data, lda, tau, &mut work, lwork, &mut info);
        }

        A
    }

    pub fn get_R(&self) -> Matrix {
        let qr = &self.fact_mat;
        let row = min(qr.row, qr.col);
        let col = qr.col;
        let mut result = matrix(vec![0f64; row * col], row, col, qr.shape);
        for i in 0..row {
            for j in i..col {
                result[(i, j)] = qr[(i, j)];
            }
        }
        result
    }
}

#[allow(non_snake_case)]
impl DPOTRF {
    pub fn get_U(&self) -> Option<Matrix> {
        if self.uplo == UPLO::Lower {
            return None;
        }

        let mat = &self.fact_mat;
        let n = mat.col;
        let mut result = matrix(vec![0f64; n.pow(2)], n, n, mat.shape);
        for i in 0 .. n {
            for j in i .. n {
                result[(i, j)] = mat[(i, j)];
            }
        }
        Some(result)
    }

    pub fn get_L(&self) -> Option<Matrix> {
        if self.uplo == UPLO::Upper {
            return None;
        }

        let mat = &self.fact_mat;
        let n = mat.col;
        let mut result = matrix(vec![0f64; n.pow(2)], n, n, mat.shape);

        for i in 0 .. n {
            for j in 0 .. i+1 {
                result[(i, j)] = mat[(i, j)];
            }
        }
        Some(result)
    }
}

#[allow(non_snake_case)]
pub fn gen_householder(a: &Vec<f64>) -> Matrix {
    let mut v = a.fmap(|t| t / (a[0] + a.norm(Norm::L2) * a[0].signum()));
    v[0] = 1f64;
    let mut H = eye(a.len());
    let vt: Matrix = v.clone().into();
    H = H - 2f64 / v.dot(&v) * (&vt * &vt.t());
    H
}

/// LU via Gaussian Elimination with Partial Pivoting
#[allow(dead_code)]
fn gepp(m: &mut Matrix) -> Vec<usize> {
    let mut r = vec![0usize; m.col - 1];
    for k in 0..(m.col - 1) {
        // Find the pivot row
        let r_k = m
            .col(k)
            .into_iter()
            .skip(k)
            .enumerate()
            .max_by(|x1, x2| x1.1.abs().partial_cmp(&x2.1.abs()).unwrap())
            .unwrap()
            .0
            + k;
        r[k] = r_k;

        // Interchange the rows r_k and k
        for j in k..m.col {
            unsafe {
                std::ptr::swap(&mut m[(k, j)], &mut m[(r_k, j)]);
                println!("Swap! k:{}, r_k:{}", k, r_k);
            }
        }
        // Form the multipliers
        for i in k + 1..m.col {
            m[(i, k)] = -m[(i, k)] / m[(k, k)];
        }
        // Update the entries
        for i in k + 1..m.col {
            for j in k + 1..m.col {
                m[(i, j)] += m[(i, k)] * m[(k, j)];
            }
        }
    }
    r
}

/// LU via Gauss Elimination with Complete Pivoting
fn gecp(m: &mut Matrix) -> (Vec<usize>, Vec<usize>) {
    let n = m.col;
    let mut r = vec![0usize; n - 1];
    let mut s = vec![0usize; n - 1];
    for k in 0..n - 1 {
        // Find pivot
        let (r_k, s_k) = match m.shape {
            Col => {
                let mut row_ics = 0usize;
                let mut col_ics = 0usize;
                let mut max_val = 0f64;
                for i in k..n {
                    let c = m
                        .col(i)
                        .into_iter()
                        .skip(k)
                        .enumerate()
                        .max_by(|x1, x2| x1.1.abs().partial_cmp(&x2.1.abs()).unwrap())
                        .unwrap();
                    let c_ics = c.0 + k;
                    let c_val = c.1.abs();
                    if c_val > max_val {
                        row_ics = c_ics;
                        col_ics = i;
                        max_val = c_val;
                    }
                }
                (row_ics, col_ics)
            }
            Row => {
                let mut row_ics = 0usize;
                let mut col_ics = 0usize;
                let mut max_val = 0f64;
                for i in k..n {
                    let c = m
                        .row(i)
                        .into_iter()
                        .skip(k)
                        .enumerate()
                        .max_by(|x1, x2| x1.1.abs().partial_cmp(&x2.1.abs()).unwrap())
                        .unwrap();
                    let c_ics = c.0 + k;
                    let c_val = c.1.abs();
                    if c_val > max_val {
                        col_ics = c_ics;
                        row_ics = i;
                        max_val = c_val;
                    }
                }
                (row_ics, col_ics)
            }
        };
        r[k] = r_k;
        s[k] = s_k;

        // Interchange rows
        for j in k..n {
            unsafe {
                std::ptr::swap(&mut m[(k, j)], &mut m[(r_k, j)]);
            }
        }

        // Interchange cols
        for i in 0..n {
            unsafe {
                std::ptr::swap(&mut m[(i, k)], &mut m[(i, s_k)]);
            }
        }

        // Form the multipliers
        for i in k + 1..n {
            m[(i, k)] = -m[(i, k)] / m[(k, k)];
            for j in k + 1..n {
                m[(i, j)] += m[(i, k)] * m[(k, j)];
            }
        }
    }
    (r, s)
}
