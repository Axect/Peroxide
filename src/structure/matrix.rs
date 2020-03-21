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
//! * Description: Same as R - `matrix(Vector, Row, Col, Shape)`
//! * Type: `matrix(Vec<T>, usize, usize, Shape) where T: std::convert::Into<f64> + Copy`
//!     * `Shape`: `Enum` for matrix shape - `Row` & `Col`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
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
//!     use peroxide::*;
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
//!     use peroxide::*;
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
//!     a = matrix(1:4:1, 2, 2, Row)
//!     print(a)
//!     #      [,1] [,2]
//!     # [1,]    1    2
//!     # [2,]    3    4
//!     ```
//!
//! * For `Peroxide`,
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
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
//!     extern crate peroxide;
//!     use peroxide::*;
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
//! ### CSV (Legacy)
//!
//! * `write(&self, file_path: &str)`: Write matrix to csv
//! * `write_with_header(&self, file_path, header: Vec<&str>)`: Write with header
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = ml_matrix("1 2;3 4");
//!         a.write("example_data/matrix.csv").expect("Can't write file");
//!
//!         let b = ml_matrix("1 2; 3 4; 5 6");
//!         b.write_with_header("example_data/header.csv", vec!["odd", "even"])
//!             .expect("Can't write header file");
//!     }
//!     ```
//!
//! Also, you can read matrix from csv.
//!
//! * Type: `read(&str, bool, char) -> Result<Matrix, Box<Error>>`
//! * Description: `read(file_path, is_header, delimiter)`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         //let a = Matrix::read("example_data/matrix.csv", false, ',')
//!         //    .expect("Can't read matrix.csv file");
//!         //a.print();
//!         ////       c[0] c[1]
//!         //// r[0]     1    2
//!         //// r[1]     3    4
//!     }
//!     ```
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
//!     use peroxide::*;
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
//!     extern crate peroxide;
//!     use peroxide::*;
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
//!     use peroxide::*;
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
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let mut a = matrix!(1;4;1, 2, 2, Row);
//!         a[(0,0)].print(); // 1
//!         a[(0,0)] = 2f64; // Modify component
//!         a.print();
//!         //       c[0] c[1]
//!         //  r[0]    2    2
//!         //  r[1]    3    4
//!     }
//!     ```
//!
//! ## Conversion to vector
//!
//! * Just use `row` or `col` method (I already showed at Basic method section).
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         a.row(0).print(); // [1, 2]
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
//!     use peroxide::*;
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
//!     extern crate peroxide;
//!     use peroxide::*;
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
//! * `lu` returns `Option<PQLU>`
//!     * `PQLU` has four field - `p`, `q`, `l` , `u`
//!     * `p` means row permutations
//!     * `q` means column permutations
//!     * `l` means lower triangular matrix
//!     * `u` menas upper triangular matrix
//! * The structure of `PQLU` is as follows:
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     #[derive(Debug, Clone)]
//!     pub struct PQLU {
//!         pub p: Perms,
//!         pub q: Perms,
//!         pub l: Matrix,
//!         pub u: Matrix,
//!     }
//!
//!     pub type Perms = Vec<(usize, usize)>;
//!     ```
//!
//! * Example of LU decomposition:
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = matrix(c!(1,2,3,4), 2, 2, Row);
//!         let pqlu = a.lu().unwrap(); // unwrap because of Option
//!         let (p,q,l,u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);
//!         assert_eq!(p, vec![(0,1)]); // swap 0 & 1 (Row)
//!         assert_eq!(q, vec![(0,1)]); // swap 0 & 1 (Col)
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
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         assert_eq!(a.det(), -2f64);
//!     }
//!     ```
//!
//! ## Inverse matrix
//!
//! * Peroxide uses LU decomposition to obtain inverse matrix.
//! * It needs two sub functions - `inv_l`, `inv_u`
//!     * For inverse of `L, U`, I use block partitioning. For example, for lower triangular matrix :
//!     $$ \begin{aligned} L &= \begin{pmatrix} L_1 & \mathbf{0} \\\ L_2 & L_3 \end{pmatrix} \\\ L^{-1} &= \begin{pmatrix} L_1^{-1} & \mathbf{0} \\\ -L_3^{-1}L_2 L_1^{-1} & L_3^{-1} \end{pmatrix} \end{aligned} $$
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         a.inv().unwrap().print();
//!         //      c[0] c[1]
//!         // r[0]   -2    1
//!         // r[1]  1.5 -0.5
//!     }
//!     ```
//!
//! ## Moore-Penrose Pseudo Inverse
//!
//! * $ X^\dagger = \left(X^T X\right)^{-1} X $
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = matrix!(1;4;1, 2, 2, Row);
//!         let pinv_a = a.pseudo_inv().unwrap();
//!         let inv_a = a.inv().unwrap();
//!
//!         assert_eq!(inv_a, pinv_a); // Nearly equal (not actually equal)
//!     }
//!     ```

extern crate csv;

#[cfg(feature = "O3")]
extern crate blas;
#[cfg(feature = "O3")]
extern crate lapack;

#[cfg(feature = "O3")]
use blas::{daxpy, dgemm, dgemv};
#[cfg(feature = "O3")]
use lapack::{dgecon, dgeqrf, dgetrf, dgetri, dgetrs, dorgqr};
#[cfg(feature = "O3")]
use std::f64::NAN;

use self::csv::{ReaderBuilder, StringRecord, WriterBuilder};
pub use self::Norm::*;
pub use self::Shape::{Col, Row};
use std::cmp::{max, min};
use std::convert;
pub use std::error::Error;
use std::fmt;
use std::ops::{Add, Index, IndexMut, Mul, Neg, Sub, Div};
#[allow(unused_imports)]
use structure::vector::*;
use util::useful::*;
use MutMatrix;

pub type Perms = Vec<(usize, usize)>;

/// To select matrices' binding.
///
/// Row - Row binding
/// Col - Column binding
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4], 2, 2, Row);
/// let b = matrix(vec![1,2,3,4], 2, 2, Col);
/// println!("{}", a); // Similar to [[1,2],[3,4]]
/// println!("{}", b); // Similar to [[1,3],[2,4]]
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Shape {
    Row,
    Col,
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

/// R-like matrix structure
///
/// # Examples
///
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = Matrix {
///     data: vec![1f64,2f64,3f64,4f64],
///     row: 2,
///     col: 2,
///     shape: Row,
/// }; // [[1,2],[3,4]]
/// ```
#[derive(Debug, Clone)]
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(c!(1,2,3,4), 2, 2, Row);
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = py_matrix(vec![c!(1,2), c!(3,4)]);
/// let b = matrix(c!(1,2,3,4), 2, 2, Row);
/// assert_eq!(a, b);
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = ml_matrix("1 2; 3 4");
/// let b = matrix(c!(1,2,3,4), 2, 2, Row);
/// assert_eq!(a, b);
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
    pub fn ptr(&self) -> *const f64 {
        &self.data[0] as *const f64
    }

    pub fn mut_ptr(&mut self) -> *mut f64 {
        &mut self.data[0] as *mut f64
    }

    /// Change Bindings
    ///
    /// `Row` -> `Col` or `Col` -> `Row`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
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
                matrix(data,r,c,Col)
            }
            Col => {
                for i in 0..l {
                    let s = (i * r) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                matrix(data,r,c,Row)
            }
        }
    }

    /// Spread data(1D vector) to 2D formatted String
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
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
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(c!(1,2,3,4), 2, 2, Row);
    /// assert_eq!(a.col(0), c!(1,3));
    /// ```
    pub fn col(&self, index: usize) -> Vector {
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
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(c!(1,2,3,4), 2, 2, Row);
    /// assert_eq!(a.row(0), c!(1,2));
    /// ```
    pub fn row(&self, index: usize) -> Vector {
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
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.diag(), c!(1,4));
    /// ```
    pub fn diag(&self) -> Vector {
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

    /// Write to CSV
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// a.write("example_data/test.csv");
    /// ```
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
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// a.write_round("example_data/test.csv", 0);
    /// ```
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
    /// extern crate peroxide;
    /// use peroxide::*;
    /// use std::process;
    ///
    /// let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// a.write_round("example_data/test.csv", 0);
    ///
    /// let b = Matrix::read("example_data/test.csv", false, ','); // header = false, delimiter = ','
    /// match b {
    ///     Ok(mat) => println!("{}", mat),
    ///     Err(err) => {
    ///         println!("{}", err);
    ///         process::exit(1);
    ///     }
    /// }
    /// ```
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
}

// =============================================================================
// Common Properties of Matrix & Vector
// =============================================================================

/// Common trait for Matrix & Vector
pub trait LinearOps {
    type Operator;
    fn from_matrix(m: Matrix) -> Self;
    fn to_matrix(&self) -> Matrix;
    fn transpose(&self) -> Self::Operator;
    fn t(&self) -> Self::Operator;
}

impl LinearOps for Matrix {
    type Operator = Self;

    /// Just Clone
    fn from_matrix(m: Matrix) -> Self {
        m
    }

    /// Just clone
    fn to_matrix(&self) -> Self {
        self.clone()
    }

    /// Transpose
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(vec![1,2,3,4], 2, 2, Row);
    /// println!("{}", a); // [[1,3],[2,4]]
    /// ```
    fn transpose(&self) -> Self {
        match self.shape {
            Row => matrix(self.data.clone(), self.col, self.row, Col),
            Col => matrix(self.data.clone(), self.col, self.row, Row),
        }
    }

    /// R-like transpose function
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.transpose(), a.t());
    /// ```
    fn t(&self) -> Self {
        self.transpose()
    }
}

impl LinearOps for Vec<f64> {
    type Operator = Matrix;

    fn from_matrix(m: Matrix) -> Self {
        m.data
    }

    /// Vector to Column Matrix
    fn to_matrix(&self) -> Matrix {
        let l = self.len();
        matrix(self.clone(), l, 1, Col)
    }

    /// Vector to Row Matrix
    fn transpose(&self) -> Matrix {
        let l = self.len();
        matrix(self.clone(), 1, l, Row)
    }

    /// R-like Syntax
    fn t(&self) -> Matrix {
        self.transpose()
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

/// Element-wise addition of &Matrix
impl<'a, 'b> Add<&'b Matrix> for &'a Matrix {
    type Output = Matrix;

    fn add(self, other: &'b Matrix) -> Self::Output {
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
}

/// Element-wise addition between Matrix & f64
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// assert_eq!(a + 1, matrix!(2;5;1, 2, 2, Row));
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// assert_eq!(1f64 + a, matrix!(2;5;1, 2, 2, Row));
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// assert_eq!(1 + a, matrix!(2;5;1, 2, 2, Row));
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// assert_eq!(1 as usize + a, matrix!(2;5;1, 2, 2, Row));
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
/// use peroxide::*;
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
/// use peroxide::*;
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

/// Subtraction between &'a Matrix
impl<'a, 'b> Sub<&'b Matrix> for &'a Matrix {
    type Output = Matrix;

    fn sub(self, other: &'b Matrix) -> Self::Output {
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4],2,2,Row);
/// assert_eq!(a - 1f64, matrix!(0;3;1, 2, 2, Row));
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

impl<'a> Mul<f64> for &'a Matrix {
    type Output = Matrix;

    fn mul(self, other: f64) -> Self::Output {
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

impl<'a> Mul<i64> for &'a Matrix {
    type Output = Matrix;

    fn mul(self, other: i64) -> Self::Output {
        self.mul(other as f64)
    }
}

impl<'a> Mul<i32> for &'a Matrix {
    type Output = Matrix;

    fn mul(self, other: i32) -> Self::Output {
        self.mul(other as f64)
    }
}

impl<'a> Mul<usize> for &'a Matrix {
    type Output = Matrix;

    fn mul(self, other: usize) -> Self::Output {
        self.mul(other as f64)
    }
}

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
        other.mul(self)
    }
}

impl<'a> Mul<&'a Matrix> for i64 {
    type Output = Matrix;

    fn mul(self, other: &'a Matrix) -> Matrix {
        other.mul(self as f64)
    }
}

impl<'a> Mul<&'a Matrix> for i32 {
    type Output = Matrix;

    fn mul(self, other: &'a Matrix) -> Matrix {
        other.mul(self as f64)
    }
}

impl<'a> Mul<&'a Matrix> for usize {
    type Output = Matrix;

    fn mul(self, other: &'a Matrix) -> Matrix {
        other.mul(self as f64)
    }
}

/// Matrix Multiplication
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// let b = matrix!(1;4;1, 2, 2, Col);
/// assert_eq!(a * b, matrix(c!(5, 11, 11, 25), 2, 2, Row));
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
impl Mul<Vector> for Matrix {
    type Output = Self;

    fn mul(self, other: Vector) -> Self::Output {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let x = &other;
                let mut y = vec![0f64; self.row];
                let A = &self.data;
                let m_i32 = self.row as i32;
                let n_i32 = self.col as i32;
                match self.shape {
                    Row => unsafe {
                        dgemv(b'T', m_i32, n_i32, 1f64, A, n_i32, x, 1, 0f64, &mut y, 1);
                    },
                    Col => unsafe {
                        dgemv(b'N', m_i32, n_i32, 1f64, A, m_i32, x, 1, 0f64, &mut y, 1);
                    },
                }
                matrix(y, self.row, 1, self.shape)
            }
            _ => self.mul(other.to_matrix()),
        }
    }
}

#[allow(non_snake_case)]
impl<'a, 'b> Mul<&'b Vector> for &'a Matrix {
    type Output = Matrix;

    fn mul(self, other: &'b Vector) -> Self::Output {
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
                        dgemv(b'T', m_i32, n_i32, 1f64, A, n_i32, x, 1, 0f64, &mut y, 1);
                    },
                    Col => unsafe {
                        dgemv(b'N', m_i32, n_i32, 1f64, A, m_i32, x, 1, 0f64, &mut y, 1);
                    },
                }
                matrix(y, self.row, 1, self.shape)
            }
            _ => self.mul(&other.to_matrix()),
        }
    }
}

/// Matrix multiplication for Vector vs Matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// let v = c!(1,2);
/// assert_eq!(v * a, matrix(c!(7,10),1,2,Row));
/// ```
impl Mul<Matrix> for Vector {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Self::Output {
        assert_eq!(self.len(), other.row);
        let l = self.len();
        matrix(self, 1, l, Col).mul(other)
    }
}

impl<'a, 'b> Mul<&'b Matrix> for &'a Vector {
    type Output = Matrix;

    fn mul(self, other: &'b Matrix) -> Self::Output {
        assert_eq!(self.len(), other.row);
        let l = self.len();
        matrix(self.clone(), 1, l, Col).mul(other.clone())
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
                    daxpy(n_i32, 1f64/a_f64, x, 1, &mut y, 1);
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
/// use peroxide::*;
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
            Row => {
                unsafe {
                    &*p.add(i * self.col + j)
                }
            }   
            Col => {
                unsafe {
                    &*p.add(i + j * self.row)
                }  
            },
        }
    }
}

/// IndexMut for Matrix (Assign)
///
/// `(usize, usize) -> f64`
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let mut a = matrix!(1;4;1, 2, 2, Row);
/// a[(1,1)] = 10.0;
/// assert_eq!(a, matrix(c!(1,2,3,10), 2, 2, Row));
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
                unsafe {
                    &mut *p.add(idx)
                }
            }
            Col => {
                let idx = i + j * r;
                unsafe {
                    &mut *p.add(idx)
                }
            }
        }
    }
}

// =============================================================================
// Functional Programming Tools (Hand-written)
// =============================================================================
pub trait FP {
    fn take_row(&self, n: usize) -> Matrix;
    fn take_col(&self, n: usize) -> Matrix;
    fn skip_row(&self, n: usize) -> Matrix;
    fn skip_col(&self, n: usize) -> Matrix;
    fn fmap<F>(&self, f: F) -> Matrix
    where
        F: Fn(f64) -> f64;
    fn col_map<F>(&self, f: F) -> Matrix
    where
        F: Fn(Vec<f64>) -> Vec<f64>;
    fn row_map<F>(&self, f: F) -> Matrix
    where
        F: Fn(Vec<f64>) -> Vec<f64>;
    fn col_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<f64>) -> Vec<f64>;
    fn row_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<f64>) -> Vec<f64>;
    fn reduce<F, T>(&self, init: T, f: F) -> f64
    where
        F: Fn(f64, f64) -> f64,
        T: convert::Into<f64>;
    fn zip_with<F>(&self, f: F, other: &Matrix) -> Matrix
    where
        F: Fn(f64, f64) -> f64;
}

impl FP for Matrix {
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
        let result = self
            .data
            .clone()
            .into_iter()
            .map(|x| f(x))
            .collect::<Vec<f64>>();
        matrix(result, self.row, self.col, self.shape)
    }

    fn col_map<F>(&self, f: F) -> Matrix
    where
        F: Fn(Vec<f64>) -> Vec<f64>,
    {
        let mut result = matrix(vec![0f64; self.row * self.col], self.row, self.col, Col);

        for i in 0..self.row {
            result.subs_col(i, &f(self.col(i)));
        }

        result
    }

    fn row_map<F>(&self, f: F) -> Matrix
    where
        F: Fn(Vec<f64>) -> Vec<f64>,
    {
        let mut result = matrix(vec![0f64; self.row * self.col], self.row, self.col, Row);

        for i in 0..self.col {
            result.subs_row(i, &f(self.row(i)));
        }

        result
    }

    fn col_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<f64>) -> Vec<f64>,
    {
        for i in 0 .. self.col {
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
        for i in 0 .. self.col {
            unsafe {
                let mut p = self.row_mut(i);
                let fv = f(self.row(i));
                for j in 0 .. p.len() {
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
        self.data
            .clone()
            .into_iter()
            .fold(init.into(), |x, y| f(x, y))
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
            .clone()
            .into_iter()
            .zip(a.data)
            .map(|(x, y)| f(x, y))
            .collect::<Vec<f64>>();
        matrix(result, self.row, self.col, self.shape)
    }
}

// =============================================================================
// Linear Algebra
// =============================================================================

/// Norm Enum
#[derive(Debug, Copy, Clone)]
pub enum Norm {
    Frobenius,
    PQ(usize, usize),
    One,
    Infinity,
}

/// Linear algebra trait
pub trait LinearAlgebra {
    fn norm(&self, norm: Norm) -> f64;
    fn lu(&self) -> Option<PQLU>;
    fn det(&self) -> f64;
    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix);
    fn inv(&self) -> Option<Matrix>;
    fn pseudo_inv(&self) -> Option<Matrix>;
}

pub fn diag(n: usize) -> Matrix {
    let mut v: Vec<f64> = vec![0f64; n * n];
    for i in 0..n {
        let idx = i * (n + 1);
        v[idx] = 1f64;
    }
    matrix(v, n, n, Row)
}

/// Data structure for LU decomposition
///
/// # Usage
/// ```rust
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = ml_matrix("1 2;3 4");
/// let pqlu = a.lu().unwrap();
/// let (p, q, l, u) = pqlu.extract();
/// // p, q are permutations
/// // l, u are matrices
/// l.print(); // lower triangular
/// u.print(); // upper triangular
/// ```
#[derive(Debug, Clone)]
pub struct PQLU {
    pub p: Perms,
    pub q: Perms,
    pub l: Matrix,
    pub u: Matrix,
}

impl PQLU {
    pub fn extract(&self) -> (Perms, Perms, Matrix, Matrix) {
        (
            self.p.clone(),
            self.q.clone(),
            self.l.clone(),
            self.u.clone(),
        )
    }
}

impl LinearAlgebra for Matrix {
    /// Matrix norm
    ///
    /// # Kinds
    /// * `Frobenius` : Frobenius norm
    /// * `PQ(usize, usize)` : L_pq norm
    /// * `One` : 1-norm
    /// * `Infinity` : Infinity norm
    fn norm(&self, norm: Norm) -> f64 {
        match norm {
            Frobenius => {
                let mut s = 0f64;
                for i in 0..self.data.len() {
                    s += self.data[i].powi(2);
                }
                s.sqrt()
            }
            PQ(p, q) => {
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
            One => {
                let mut m = std::f64::MIN;
                match self.shape {
                    Row => self.change_shape().norm(One),
                    Col => {
                        for c in 0..self.col {
                            let s = self.col(c).reduce(0f64, |x, y| x + y);
                            if s > m {
                                m = s;
                            }
                        }
                        m
                    }
                }
            }
            Infinity => {
                let mut m = std::f64::MIN;
                let a = match self.shape {
                    Row => self.clone(),
                    Col => self.change_shape(),
                };
                for r in 0..a.row {
                    let s = a.row(r).reduce(0f64, |x, y| x + y);
                    if s > m {
                        m = s;
                    }
                }
                m
            }
        }
    }

    /// LU Decomposition Implements
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
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(vec![1,2,3,4], 2, 2, Row);
    /// let pqlu = a.lu().unwrap();
    /// let (p,q,l,u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);
    /// assert_eq!(p, vec![(0,1)]); // swap 0 & 1 (Row)
    /// assert_eq!(q, vec![(0,1)]); // swap 0 & 1 (Col)
    /// assert_eq!(l, matrix(c!(1,0,0.5,1),2,2,Row));
    /// assert_eq!(u, matrix(c!(4,3,0,-0.5),2,2,Row));
    /// ```
    fn lu(&self) -> Option<PQLU> {
        assert_eq!(self.col, self.row);
        let n = self.row;
        let len: usize = n * n;

        let mut l: Self = matrix(vec![0f64; len], n, n, self.shape);
        let mut u: Self = matrix(vec![0f64; len], n, n, self.shape);

        // ---------------------------------------
        // Pivoting - Complete
        // ---------------------------------------
        // Permutations
        let mut p: Perms = Vec::new();
        let mut q: Perms = Vec::new();

        let mut container = self.clone();

        for k in 0..(n - 1) {
            // Initialize maximum & Position
            let mut m = 0f64;
            let mut row_idx: usize = k;
            let mut col_idx: usize = k;
            // Find Maximum value & Position
            for i in k..n {
                for j in k..n {
                    let temp = container[(i, j)];
                    if temp.abs() > m.abs() {
                        m = temp;
                        row_idx = i;
                        col_idx = j;
                    }
                }
            }
            if k != row_idx {
                unsafe {
                    container.swap(k, row_idx, Row); // Row perm
                }
                p.push((k, row_idx));
            }
            if k != col_idx {
                unsafe {
                    container.swap(k, col_idx, Col); // Col perm
                }
                q.push((k, col_idx));
            }
        }

        // ---------------------------------------
        // Obtain L & U
        // ---------------------------------------
        let reference = container.clone();

        // Initialize U
        for i in 0..n {
            u[(0, i)] = reference[(0, i)];
        }

        // Initialize L
        for i in 0..n {
            l[(i, i)] = 1f64;
        }

        for i in 0..n {
            for k in i..n {
                let mut s = 0f64;
                for j in 0..i {
                    s += l[(i, j)] * u[(j, k)];
                }
                u[(i, k)] = reference[(i, k)] - s;
                // Check non-zero diagonal
                if nearly_eq(u[(i, i)], 0) {
                    return None;
                }
            }

            for k in (i + 1)..n {
                let mut s = 0f64;
                for j in 0..i {
                    s += l[(k, j)] * u[(j, i)];
                }
                l[(k, i)] = (reference[(k, i)] - s) / u[(i, i)];
            }
        }

        Some(PQLU { p, q, l, u })
    }

    /// Determinant
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.det(), -2f64);
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
            _ => {
                match self.lu() {
                    None => 0f64,
                    Some(pqlu) => {
                        let (p, q, _l, u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);

                        // sgn of perms
                        let sgn_p = 2.0 * (p.len() % 2) as f64 - 1.0;
                        let sgn_q = 2.0 * (q.len() % 2) as f64 - 1.0;

                        u.diag().reduce(1f64, |x, y| x * y) * sgn_p * sgn_q
                    }
                }
            }
        }
    }

    /// Block Partition
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;16;1, 4, 4, Row);
    /// let (m1, m2, m3, m4) = a.block();
    /// assert_eq!(m1, matrix(c!(1,2,5,6), 2, 2, Row));
    /// assert_eq!(m2, matrix(c!(3,4,7,8), 2, 2, Row));
    /// assert_eq!(m3, matrix(c!(9,10,13,14), 2, 2, Row));
    /// assert_eq!(m4, matrix(c!(11,12,15,16), 2, 2, Row));
    ///
    /// let b = matrix!(1;16;1, 4, 4, Col);
    /// let (m1, m2, m3, m4) = b.block();
    /// assert_eq!(m1, matrix(c!(1,2,5,6), 2, 2, Col));
    /// assert_eq!(m3, matrix(c!(3,4,7,8), 2, 2, Col));
    /// assert_eq!(m2, matrix(c!(9,10,13,14), 2, 2, Col));
    /// assert_eq!(m4, matrix(c!(11,12,15,16), 2, 2, Col));
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
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// // Non-singular
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.inv().unwrap(), matrix(c!(-2,1,1.5,-0.5),2,2,Row));
    ///
    /// // Singular
    /// let b = matrix!(1;9;1, 3, 3, Row);
    /// assert_eq!(b.inv(), None);
    /// ```
    fn inv(&self) -> Option<Self> {
        match () {
            #[cfg(feature = "O3")]
            () => {
                let opt_dgrf = lapack_dgetrf(self);
                match opt_dgrf {
                    None => None,
                    Some(dgrf) => lapack_dgetri(&dgrf),
                }
            }
            _ => match self.lu() {
                None => None,
                Some(pqlu) => {
                    let (p, q, l, u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);
                    let mut m = inv_u(u) * inv_l(l);
                    for (idx1, idx2) in q.into_iter().rev() {
                        unsafe {
                            m.swap(idx1, idx2, Row);
                        }
                    }
                    for (idx1, idx2) in p.into_iter().rev() {
                        unsafe {
                            m.swap(idx1, idx2, Col);
                        }
                    }
                    Some(m)
                }
            },
        }
    }

    /// Moore-Penrose Pseudo inverse
    ///
    /// # Description
    /// `(X^T X)^{-1} X`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// let inv_a = a.inv().unwrap();
    /// let pse_a = a.pseudo_inv().unwrap();
    ///
    /// assert_eq!(inv_a, pse_a); // Nearly equal
    /// ```
    fn pseudo_inv(&self) -> Option<Self> {
        let xtx = &self.t() * self;
        let inv_temp = xtx.inv();

        match inv_temp {
            None => None,
            Some(m) => Some(m * self.t()),
        }
    }
}

#[allow(non_snake_case)]
pub fn solve(A: &Matrix, b: &Matrix) -> Option<Matrix> {
    match () {
        #[cfg(feature = "O3")]
        () => {
            let opt_dgrf = lapack_dgetrf(A);
            match opt_dgrf {
                None => None,
                Some(dgrf) => match dgrf.status {
                    LAPACK_STATUS::Singular => None,
                    LAPACK_STATUS::NonSingular => lapack_dgetrs(&dgrf, b),
                },
            }
        }
        _ => {
            let opt_a_inv = A.inv();
            match opt_a_inv {
                None => None,
                Some(a_inv) => Some(&a_inv * b),
            }
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;16;1, 4, 4, Row);
/// let (m1, m2, m3, m4) = a.block();
/// let m = combine(m1,m2,m3,m4);
/// assert_eq!(m, a);
///
/// let b = matrix!(1;16;1, 4, 4, Col);
/// let (n1, n2, n3, n4) = b.block();
/// let n = combine(n1,n2,n3,n4);
/// assert_eq!(n, b);
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(c!(1,0,2,1), 2, 2, Row);
/// assert_eq!(inv_l(a), matrix(c!(1,0,-2,1), 2, 2, Row));
///
/// let b = matrix(c!(1,0,0,2,1,0,4,3,1), 3, 3, Row);
/// assert_eq!(inv_l(b), matrix(c!(1,0,0,-2,1,0,2,-3,1), 3, 3, Row));
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
/// extern crate peroxide;
/// use peroxide::*;
///
/// let u = matrix(c!(2,2,0,1), 2, 2, Row);
/// assert_eq!(inv_u(u), matrix(c!(0.5,-1,0,1), 2, 2, Row));
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

fn matmul(a: &Matrix, b: &Matrix) -> Matrix {
    match (a.row, a.col) {
        (p, q) if p <= 100 && q <= 100 => {
            let r_self = a.row;
            let c_self = a.col;
            let new_other = b;
            let r_other = new_other.row;
            let c_other = new_other.col;

            assert_eq!(c_self, r_other);

            let r_new = r_self;
            let c_new = c_other;

            let mut result = matrix(vec![0f64; r_new * c_new], r_new, c_new, a.shape);

            for i in 0..r_new {
                for j in 0..c_new {
                    let mut s = 0f64;
                    for k in 0..c_self {
                        s += a[(i, k)] * new_other[(k, j)];
                    }
                    result[(i, j)] = s;
                }
            }
            result
        }
        _ => {
            let (a1, a2, a3, a4) = a.block();
            let (b1, b2, b3, b4) = b.block();

            let m1 = matmul(&a1, &b1) + matmul(&a2, &b3);
            let m2 = matmul(&a1, &b2) + matmul(&a2, &b4);
            let m3 = matmul(&a3, &b1) + matmul(&a4, &b3);
            let m4 = matmul(&a3, &b2) + matmul(&a4, &b4);

            combine(m1, m2, m3, m4)
        }
    }
}

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
        let anorm = A.norm(One);
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
