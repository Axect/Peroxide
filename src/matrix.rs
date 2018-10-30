use std::convert;
use std::fmt;
use std::ops::{Add, Neg, Sub, Mul, Rem, Index, IndexMut};
pub use self::Shape::{Row, Col};
pub use vector::*;
use std::f64::{MAX, MIN};
use std::cmp::{max, min};

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
/// let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
/// let b = Matrix::new(vec![1,2,3,4], 2, 2, Col);
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

/// To convert generic vector to Matrix
pub trait CreateMatrix<T: convert::Into<f64>> {
    fn new(v: Vec<T>, x:usize, y:usize, shape: Shape) -> Matrix;
}

impl<T> CreateMatrix<T> for Matrix where T: convert::Into<f64> {
    /// Matrix generic constructor
    ///
    /// You can use any numeric type vector
    /// e.g. `u32`, `i32`, `i64`, ...
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = Matrix::new(
    ///     vec![1,2,3,4], // Can use any Into<f64> type
    ///     2,
    ///     2,
    ///     Row,
    /// );
    /// ```
    fn new(v: Vec<T>, x: usize, y: usize, shape: Shape) -> Matrix {
        Matrix {
            data: v.into_iter().map(|x| x.into()).collect(),
            row: x,
            col: y,
            shape: shape,
        }
    }
}

/// R-like matrix constructor
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(c!(1,2,3,4), 2, 2, Row);
/// ```
pub fn matrix<T>(v: Vec<T>, x:usize, y:usize, shape: Shape) -> Matrix where T: convert::Into<f64> {
    Matrix::new(v, x, y, shape)
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
            self.data.clone()
                .into_iter()
                .zip(other.data.clone())
                .all(|(x, y)| nearly_eq(x,y)) && self.row == other.row
        } else {
            self.eq(&other.change_shape())
        }
    }
}

/// Main matrix structure
#[allow(dead_code)]
impl Matrix {
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
    pub fn change_shape(&self) -> Matrix {
        let r = self.row;
        let c = self.col;
        assert_eq!(r*c, self.data.len());
        let l = r * c - 1;
        let mut data: Vec<f64> = self.data.clone();
        let ref_data: Vec<f64> = self.data.clone();

        match self.shape {
            Row => {
                for i in 0 .. l {
                    let s = (i * c) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                Matrix {
                    data: data.clone(),
                    row: r,
                    col: c,
                    shape: Col,
                }
            },
            Col => {
                for i in 0 .. l {
                    let s = (i * r) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                Matrix {
                    data: data.clone(),
                    row: r,
                    col: c,
                    shape: Row,
                }
            },
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

        // Find maximum length of data
        let sample = self.data.clone();
        let mut space: usize = sample.into_iter()
            .map(|x|
                min(
                    format!("{:.4}", x).len(),
                    x.to_string().len(),
                ) // Choose minimum of approx vs normal
            )
            .fold(0, |x, y| max(x,y)) + 1;

        if space < 5 {
            space = 5;
        }

        let mut result = String::new();

        result.push_str(&tab("", 5));
        for i in 0 .. c {
            result.push_str(&tab(&format!("c[{}]", i), space)); // Header
        }
        result.push('\n');

        for i in 0 .. r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for j in 0 .. c {
                let st1 = format!("{:.4}",self[(i, j)]); // Round at fourth position
                let st2 = self[(i,j)].to_string();       // Normal string
                let mut st = st2.clone();

                // Select more small thing
                if st1.len() < st2.len() {
                    st = st1;
                }

                result.push_str(&tab(&st, space));
            }
            if i == (r-1) {
                break;
            }
            result.push('\n');
        }

        return result;
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
    pub fn transpose(&self) -> Matrix {
        match self.shape {
            Row => matrix(self.data.clone(), self.col, self.row, Col),
            Col => matrix(self.data.clone(), self.col, self.row, Row)
        }
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
        let mut container: Vec<f64> = Vec::new();
        match self.shape {
            Row => {
                let l: usize = self.row * self.col;
                for i in 0 .. l {
                    if i % self.col == index {
                        container.push(self.data[i]);
                    }
                }
            },
            Col => {
                let s: usize = self.row * index;
                container = self.data.clone().into_iter()
                    .skip(s)
                    .take(self.row).collect::<Vec<f64>>();
            }
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
        let mut container: Vec<f64> = Vec::new();
        match self.shape {
            Row => {
                let s: usize = self.col * index;
                container = self.data.clone().into_iter()
                    .skip(s)
                    .take(self.col).collect::<Vec<f64>>();
            },
            Col => {
                let l: usize = self.row * self.col;
                for i in 0 .. l {
                    if i % self.row == index {
                        container.push(self.data[i]);
                    }
                }
            }
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
        let mut container: Vector = Vec::new();
        let r = self.row;
        let c = self.col;
        assert_eq!(r, c);
        for i in 0 .. r {
            container.push(self.data[i * (c + 1)]);
        }
        container
    }

    /// Swap row or col
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.swap(0,1,Row), matrix(c!(3,4,1,2),2,2,Row));
    /// assert_eq!(a.swap(0,1,Col), matrix(c!(2,4,1,3),2,2,Col));
    /// ```
    pub fn swap(&self, idx1: usize, idx2: usize, shape: Shape) -> Matrix {
        match shape {
            Row => {
                let mut v: Vector = Vec::new();
                for k in 0 .. self.row {
                    if k == idx1 {
                        v.extend(&self.row(idx2));
                    } else if k == idx2 {
                        v.extend(&self.row(idx1));
                    } else {
                        v.extend(&self.row(k));
                    }
                }
                matrix(v, self.row, self.col, Row)
            },
            Col => {
                let mut v: Vector = Vec::new();
                for k in 0 .. self.col {
                    if k == idx1 {
                        v.extend(&self.col(idx2));
                    } else if k == idx2 {
                        v.extend(&self.col(idx1));
                    } else {
                        v.extend(&self.col(k));
                    }
                }
                matrix(v, self.row, self.col, Col)
            }
        }
    }
}

// =============================================================================
// Standard Operation for Matrix
// =============================================================================

/// Element-wise addition of Matrix
///
/// # Caution
/// > You should remember ownership.
/// > If you use Matrix `a,b` then you can't use them after.
impl Add<Matrix> for Matrix {
    type Output = Matrix;

    fn add(self, other: Matrix) -> Matrix {
        assert_eq!(&self.row, &other.row);
        assert_eq!(&self.col, &other.col);
        if self.shape == other.shape {
            matrix(
                self.data.clone().into_iter().zip(&other.data).map(|(x,y)| x + y).collect::<Vec<f64>>(),
                self.row,
                self.col,
                self.shape,
            )
        } else {
            self.change_shape().add(other)
        }
    }
}

/// Element-wise addition between Matrix & f64
impl Add<f64> for Matrix {
    type Output = Matrix;

    fn add(self, other: f64) -> Matrix {
        self.fmap(|x| x + other)
    }
}

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
    type Output = Matrix;
    
    fn neg(self) -> Matrix {
        matrix(
            self.data.into_iter().map(|x:f64| -x).collect::<Vec<f64>>(),
            self.row,
            self.col,
            self.shape,
        )
    }
}

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
    type Output = Matrix;

     fn sub(self, other: Matrix) -> Matrix {
        self.add(other.neg())
    }
}

impl Sub<f64> for Matrix {
    type Output = Matrix;

    fn sub(self, other: f64) -> Matrix {
        self.fmap(|x| x - other)
    }
}

/// Element-wise matrix multiplication
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4], 2, 2, Row);
/// let b = matrix(vec![1,2,3,4], 2, 2, Col);
/// println!("{}", a * b); // [[1,6],[6,16]]
/// ```
impl Mul<Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        assert_eq!(self.row, other.row);
        assert_eq!(self.col, other.col);
        self.zip_with(|x,y| x * y, &other)
    }
}

impl Mul<f64> for Matrix {
    type Output = Matrix;

    fn mul(self, other: f64) -> Matrix {
        self.fmap(|x| x * other)
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
/// assert_eq!(a % b, matrix(c!(5, 11, 11, 25), 2, 2, Row));
/// ```
impl Rem<Matrix> for Matrix {
    type Output = Matrix;

    fn rem(self, other: Matrix) -> Matrix {
        let r_self = self.row;
        let c_self = self.col;
        let r_other = other.row;
        let c_other = other.col;

        assert_eq!(c_self, r_other);

        let r_new = r_self;
        let c_new = c_other;

        let mut result = matrix(vec![0f64; r_new * c_new], r_new, c_new, self.shape);

        for i in 0 .. r_new {
            for j in 0 .. c_new {
                let mut s = 0f64;
                for k in 0 .. c_self {
                    s += self[(i, k)] * other[(k, j)];
                }
                result[(i, j)] = s;
            }
        }
        result
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
        let i = pair.0;
        let j = pair.1;
        match self.shape {
            Row => &self.data[i * self.col + j],
            Col => &self.data[i + j * self.row]
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
    fn index_mut<'a>(&'a mut self, pair: (usize, usize)) -> &'a mut f64 {
        let i = pair.0;
        let j = pair.1;
        let r = self.row;
        let c = self.col;
        match self.shape {
            Row => {
                let idx = i * c + j;
                &mut self.data[idx]
            },
            Col => {
                let idx = i + j * r;
                &mut self.data[idx]
            }
        }
    }
}

// =============================================================================
// Functional Programming Tools (Hand-written)
// =============================================================================
pub trait FP {
    fn fmap<F>(&self, f: F) -> Matrix where F: Fn(f64) -> f64;
    fn reduce<F, T>(&self, init: T, f: F) -> f64 
        where F: Fn(f64, f64) -> f64,
              T: convert::Into<f64>;
    fn zip_with<F>(&self, f: F, other: &Matrix) -> Matrix 
        where F: Fn(f64, f64) -> f64;
}

impl FP for Matrix {
    fn fmap<F>(&self, f: F) -> Matrix where F: Fn(f64) -> f64 {
        let result = self.data.clone().into_iter().map(|x| f(x)).collect::<Vec<f64>>();
        matrix(
            result,
            self.row,
            self.col,
            self.shape,
        )
    }

    fn reduce<F, T>(&self, init: T, f: F) -> f64 
        where F: Fn(f64, f64) -> f64,
              T: convert::Into<f64> {
        self.data.clone().into_iter().fold(
            init.into(),
            |x,y| f(x,y),
        )
    }

    fn zip_with<F>(&self, f: F, other: &Matrix) -> Matrix 
        where F: Fn(f64, f64) -> f64 {
        assert_eq!(self.data.len(), other.data.len());
        let mut a = other.clone();
        if self.shape != other.shape {
            a = a.change_shape();
        }
        let result = self.data.clone().into_iter()
            .zip(a.data)
            .map(|(x,y)| f(x,y))
            .collect::<Vec<f64>>();
        matrix(
            result,
            self.row,
            self.col,
            self.shape,
        )
    }
}

// =============================================================================
// Linear Algebra
// =============================================================================
/// Linear algebra trait
pub trait LinearAlgebra {
    fn lu(&self) -> Option<PQLU>;
    fn det(&self) -> f64;
    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix);
    fn inv(&self) -> Option<Matrix>;
}

pub fn diag(n: usize) -> Matrix {
    let mut v: Vec<f64> = vec![0f64; n*n];
    for i in 0 .. n {
        let idx = i * (n + 1);
        v[idx] = 1f64;
    }
    matrix(v, n, n, Row)
}

#[derive(Debug, Clone)]
pub struct PQLU {
    pub p: Perms,
    pub q: Perms,
    pub l: Matrix,
    pub u: Matrix,
}

impl LinearAlgebra for Matrix {
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

        let mut l: Matrix = matrix(vec![0f64; len], n, n, self.shape);
        let mut u: Matrix = matrix(vec![0f64; len], n, n, self.shape);

        // ---------------------------------------
        // Pivoting - Complete
        // ---------------------------------------
        // Permutations
        let mut p: Perms = Vec::new();
        let mut q: Perms = Vec::new();

        let mut container = self.clone();

        for k in 0 .. (n-1) {
            // Initialize maximum & Position
            let mut m = MIN;
            let mut row_idx: usize = k;
            let mut col_idx: usize = k;
            // Find Maximum value & Position
            for i in k .. n {
                for j in k .. n {
                    let temp = container[(i,j)];
                    if temp > m {
                        m = temp;
                        row_idx = i;
                        col_idx = j;
                    }
                }
            }
            if k != row_idx {
                container = container.swap(k, row_idx, Row); // Row perm
                p.push((k, row_idx));
            }
            if k != col_idx {
                container = container.swap(k, col_idx, Col); // Col perm
                q.push((k, col_idx));
            }
        }

        // ---------------------------------------
        // Obtain L & U
        // ---------------------------------------
        let reference = container.clone();

        // Initialize U
        for i in 0 .. n {
            u[(0, i)] = reference[(i, 0)];
        }

        // Initialize L
        for i in 0 .. n {
            l[(i, i)] = 1f64;
        }

        for i in 0 .. n {
            for k in i .. n {
                let mut s = 0f64;
                for j in 0 .. i {
                    s += l[(i, j)] * u[(j, k)];
                }
                u[(i, k)] = reference[(i, k)] - s;
                // Check non-zero diagonal
                if nearly_eq(u[(i, i)], 0) {
                    return None;
                }
            }

            for k in (i+1) .. n {
                let mut s = 0f64;
                for j in 0 .. i {
                    s += l[(k, j)] * u[(j, i)];
                }
                l[(k, i)] = (reference[(k, i)] - s) / u[(i, i)];
            }
        }

        Some(
            PQLU {
                p: p,
                q: q,
                l: l,
                u: u,
            }
        )
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
        match self.lu() {
            None => 0f64,
            Some(pqlu) => {
                let (p, q, _l, u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);

                // sgn of perms
                let sgn_p = 2.0 * (p.len() % 2) as f64 - 1.0;
                let sgn_q = 2.0 * (q.len() % 2) as f64 - 1.0;

                u.diag().reduce(1f64, |x,y| x * y) * sgn_p * sgn_q
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
    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix) {
        let r = self.row;
        let c = self.col;
        let l = min(self.row / 2, self.col / 2);
        let r_l = r - l;
        let c_l = c - l;

        let mut m1 = matrix(vec![0f64; l * l], l, l, self.shape);
        let mut m2 = matrix(vec![0f64; l * c_l], l, c_l, self.shape);
        let mut m3 = matrix(vec![0f64; r_l * l], r_l, l, self.shape);
        let mut m4 = matrix(vec![0f64; r_l * c_l], r_l, c_l, self.shape);

        for idx_row in 0 .. r {
            for idx_col in 0..c {
                match (idx_row, idx_col) {
                    (i, j) if (i < l) && (j < l) => {
                        m1[(i, j)] = self[(i, j)];
                    },
                    (i, j) if (i < l) && (j >= l) => {
                        m2[(i, j - l)] = self[(i, j)];
                    },
                    (i, j) if (i >= l) && (j < l) => {
                        m3[(i - l, j)] = self[(i, j)];
                    },
                    (i, j) if (i >= l) && (j >= l) => {
                        m4[(i - l, j - l)] = self[(i, j)];
                    },
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
    fn inv(&self) -> Option<Matrix> {
        match self.lu() {
            None => None,
            Some(pqlu) => {
                let (p, q, l, u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);
                let mut m = inv_u(u) % inv_l(l);
                for (idx1, idx2) in q.into_iter() {
                    m = m.swap(idx1, idx2, Row);
                }
                for (idx1, idx2) in p.into_iter() {
                    m = m.swap(idx1, idx2, Col);
                }
                Some(m)
            }
        }
    }
}

// =============================================================================
// Back-end Utils
// =============================================================================

#[allow(unused_comparisons)]
pub fn tab(s: &str, space: usize) -> String {
    let l = s.len();
    let mut m: String = String::new();
    let fs = format!("{}{}", " ".repeat(space - l), s);
    m.push_str(&fs);
    return m;
}

pub fn quot_rem(x: usize, y: usize) -> (i32, i32) {
    ((x / y) as i32, (x % y) as i32)
}

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
    let l = m1.col;
    let c_l = m2.col;
    let r_l = m3.row;

    let r = l + r_l;
    let c = l + c_l;

    let mut m = matrix(vec![0f64; r*c], r, c, m1.shape);

    for idx_row in 0 .. r {
        for idx_col in 0..c {
            match (idx_row, idx_col) {
                (i, j) if (i < l) && (j < l) => {
                    m[(i, j)] = m1[(i, j)];
                },
                (i, j) if (i < l) && (j >= l) => {
                    m[(i, j)] = m2[(i, j - l)];
                },
                (i, j) if (i >= l) && (j < l) => {
                    m[(i, j)] = m3[(i - l, j)];
                },
                (i, j) if (i >= l) && (j >= l) => {
                    m[(i, j)] = m4[(i - l, j - l)];
                },
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
            m[(1, 0)] = - m[(1, 0)];
            m
        },
        _ => {
            let (l1, l2, l3, l4) = l.block();

            let m1 = inv_l(l1);
            let m2 = l2;
            let m4 = inv_l(l4);
            let m3 = -((m4.clone() % l3) % m1.clone());

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
        },
        2 => {
            let a = w[(0, 0)];
            let b = w[(0, 1)];
            let c = w[(1, 1)];
            let d = a * c;

            w[(0, 0)] = 1f64 / a;
            w[(0, 1)] = - b / d;
            w[(1, 1)] = 1f64 / c;
            w
        },
        _ => {
            let (u1, u2, u3, u4) = u.block();
            let m1 = inv_u(u1);
            let m3 = u3;
            let m4 = inv_u(u4);
            let m2 = -(m1.clone() % u2 % m4.clone());

            combine(m1, m2, m3, m4)
        }
    }
}

/// near equal
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// assert!(nearly_eq(1.0/3.0 * 3.0, 1));
/// ```
pub fn nearly_eq<S, T>(x: S, y: T) -> bool
    where S: convert::Into<f64>,
          T: convert::Into<f64> {
    let mut b: bool = false;
    let e = 1e-7;
    let p: f64 = x.into().abs();
    let q: f64 = y.into().abs();
    if (p - q).abs() < e {
        b = true;
    } else if (p - q).abs() / (p + q).min(MAX) < e {
        b = true;
    }
    b
}
