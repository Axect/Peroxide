use std::convert;
use std::fmt;
use std::ops::{Add, Neg, Sub, Mul, Rem, Index};
pub use self::Shape::{Row, Col};
pub use vector_macro::*;
use std::f64::{MAX};

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
/// let a = matrix(vec![1,2,3,4], 2, 2, Row);
/// ```
pub fn matrix<T>(v: Vec<T>, x:usize, y:usize, shape: Shape) -> Matrix where T: convert::Into<f64> {
    Matrix::new(v, x, y, shape)
}

/// More R like Matrix constructor (Macro)
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row); // start;end;step
/// let b = matrix(c![1,2,3,4], 2, 2, Row);
/// let c = matrix(vec![1,2,3,4], 2, 2, Row); // Normal function
/// assert!(a == b && b == c);
/// ```
#[macro_export]
macro_rules! matrix {
    ( $start:expr;$end:expr;$step:expr,$row:expr,$col:expr,$shape:expr ) => {
        {
            matrix(
                seq!($start,$end,$step),
                $row,
                $col,
                $shape
            )
        }
    };
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
    /// let a = Matrix::new(vec![1,2,3,4],2,2,Row);
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
    /// let a = Matrix::new(vec![1,2,3,4],2,2,Row);
    /// println!("{}", a.spread()); // same as println!("{}", a);
    /// // Result:
    /// //       c[0] c[1]
    /// // r[0]     1    3
    /// // r[1]     2    4
    /// ```
    pub fn spread(&self) -> String {
        assert_eq!(self.row * self.col, self.data.len());
        let _rows = self.row;
        let cols = self.col;
        
        let mut result = String::new();
        result += &tab("");

        // Make header
        for i in 0 .. cols {
            let s = format!("c[{}]", i);
            result += &tab(&s);
            result += &s;
        }
        result += "\n";

        match self.shape {
            Row => {
                let data = self.data.clone();
                let temp: Vec<String> = data.into_iter().map(|x| x.to_string()).collect();
                let ts: Vec<String> = temp.clone().into_iter().take(cols).collect();
                let mut ss: Vec<String> = temp.into_iter().skip(cols).collect();
                let mut n: usize = 0;
                let s = format!("r[{}]", n);
                result += &s;
                result += &tab(&s);
                for txt in ts.into_iter() {
                    result += &tab(&txt);
                    result += &txt;
                }
                while ss.len() >= cols {
                    result += "\n";
                    let ts: Vec<String> = ss.clone().into_iter().take(cols).collect();
                    n += 1;
                    let s = format!("r[{}]", n);
                    result += &s;
                    result += &tab(&s);
                    let tl = ts.len();
                    for i in 0 .. tl {
                        result += &tab(&ts[i]);
                        result += &ts[i];
                    }
                    ss = ss.into_iter().skip(cols).collect();
                }
            },
            Col => {
                let mat = self.change_shape();
                let data = mat.data.clone();
                let temp: Vec<String> = data.into_iter().map(|x| x.to_string()).collect();
                let ts: Vec<String> = temp.clone().into_iter().take(cols).collect();
                let mut ss: Vec<String> = temp.into_iter().skip(cols).collect();
                let mut n: usize = 0;
                let s = format!("r[{}]", n);
                result += &s;
                result += &tab(&s);
                for txt in ts.into_iter() {
                    result += &tab(&txt);
                    result += &txt;
                }
                while ss.len() >= cols {
                    result += "\n";
                    let ts: Vec<String> = ss.clone().into_iter().take(cols).collect();
                    n += 1;
                    let s = format!("r[{}]", n);
                    result += &s;
                    result += &tab(&s);
                    let tl = ts.len();
                    for i in 0 .. tl {
                        result += &tab(&ts[i]);
                        result += &ts[i];
                    }
                    ss = ss.into_iter().skip(cols).collect();
                }
            }
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
    /// let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
    /// println!("{}", a); // [[1,3],[2,4]]
    /// ```
    pub fn transpose(&self) -> Matrix {
        match self.shape {
            Row => Matrix::new(self.data.clone(), self.col, self.row, Col),
            Col => Matrix::new(self.data.clone(), self.col, self.row, Row)
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
                    .take(self.col).collect::<Vec<f64>>();
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
                    .take(self.row).collect::<Vec<f64>>();
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
    pub fn swap(&self, i: usize, j: usize, shape: Shape) -> Matrix {
        match shape {
            Row => {
                let mut v: Vector = Vec::new();
                for k in 0 .. self.row {
                    if k == i {
                        v.extend(&self.row(j));
                    } else if k == j {
                        v.extend(&self.row(i));
                    } else {
                        v.extend(&self.row(k));
                    }
                }
                matrix(v, self.row, self.col, Row)
            },
            Col => {
                let mut v: Vector = Vec::new();
                for k in 0 .. self.col {
                    if k == i {
                        v.extend(&self.col(j));
                    } else if k == j {
                        v.extend(&self.col(i));
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
            Matrix::new(
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
/// let a = Matrix::new(vec![1,2,3,4],2,2,Row);
/// println!("{}", -a); // [[-1,-2],[-3,-4]]
/// ```
impl Neg for Matrix {
    type Output = Matrix;
    
    fn neg(self) -> Matrix {
        Matrix::new(
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
/// let a = Matrix::new(vec![1,2,3,4],2,2,Row);
/// let b = Matrix::new(vec![1,2,3,4],2,2,Col);
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
/// let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
/// let b = Matrix::new(vec![1,2,3,4], 2, 2, Col);
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
/// let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
/// let b = Matrix::new(vec![1,2,3,4], 2, 2, Col);
/// println!("{}", a % b); // [[5, 11], [11, 25]]
/// ```
impl Rem<Matrix> for Matrix {
    type Output = Matrix;

    // TODO: Not use change shape!
    fn rem(self, other: Matrix) -> Matrix {
        assert_eq!(self.col, other.row);
        if self.shape == Row && other.shape == Col {
            let mut container: Vec<f64> = Vec::new();
            for i in 0 .. self.row {
                let p = self.data.clone().into_iter().skip(self.col * i).take(self.col).collect::<Vec<f64>>();
                for j in 0 .. other.col {
                    let q = other.data.clone().into_iter().skip(other.row * j).take(other.row).collect::<Vec<f64>>();
                    let s: f64 = p.clone().into_iter().zip(&q).map(|(x, y)| x * y).fold(
                        0f64,
                        |x, y| x + y,
                    );
                    container.push(s);
                }
            }
            Matrix::new(
                container,
                self.row,
                other.col,
                Row,
            )
        } else if self.shape == Col && other.shape == Row {
            self.change_shape().rem(other.change_shape())
        } else if self.shape == Col {
            self.change_shape().rem(other)
        } else {
            self.rem(other.change_shape())
        }
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
        Matrix::new(
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
        Matrix::new(
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
    fn plu(&self) -> (Matrix, Matrix, Matrix);
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

impl LinearAlgebra for Matrix {
    /// PLU Decomposition Implements
    ///
    /// # Caution
    /// It returns tuple of matrices - (Matrix, Matrix, Matrix)
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(vec![1,2,3,4], 2, 2, Row);
    /// let plu = a.plu();
    /// //assert_eq!(plu.1, matrix(vec![1,0,3,1], 2, 2, Row));
    /// //assert_eq!(plu.1, matrix(vec![1,2,0,-2], 2, 2, Row));
    /// ```
    fn plu(&self) -> (Matrix, Matrix, Matrix) {
        assert_eq!(self.col, self.row);
        let n = self.row;
        let len: usize = n * n;
        let mut l_vec: Vec<f64> = vec![0f64; len]; // Row based
        let mut u_vec: Vec<f64> = vec![0f64; len]; // Row based

        let c1 = self.col(0);
        let p = pivot(c1);
        let reference = p.clone() % self.clone();

        // Initialize U
        match reference.shape {
            Row => {
                for i in 0 .. n {
                    u_vec[i] = reference.data[i];
                }
            },
            Col => {
                for i in 0 .. n {
                    let j: usize = i * n;
                    u_vec[j] = reference.data[j];
                }
            }
        }

        // Initialize L
        for i in 0 .. n {
            let j = i * (n + 1);
            l_vec[j] = 1f64;
        }

        for i in 0 .. n {
            for k in i .. n {
                let mut s = 0f64;
                for j in 0 .. i {
                    s += l_vec[i*n + j] * u_vec[j*n + k];
                }
                u_vec[i*n + k] = reference[(i, k)] - s;
            }

            for k in (i+1) .. n {
                let mut s = 0f64;
                for j in 0 .. i {
                    s += l_vec[k*n + j] * u_vec[j*n + i];
                }
                l_vec[k*n + i] = (reference[(k, i)] - s) / u_vec[i*n + i];
            }
        }

        let l = matrix(l_vec, n, n, Row);
        let u = matrix(u_vec, n, n, Row);

        return (p, l, u);
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
        let (p, _l, u) = self.plu();
        let sgn = to_perm(p).sign();
        u.diag().reduce(1, |x,y| x * y) * sgn
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
    /// assert_eq!(m2, matrix(c!(3,4,7,8), 2, 2, Col));
    /// assert_eq!(m3, matrix(c!(9,10,13,14), 2, 2, Col));
    /// assert_eq!(m4, matrix(c!(11,12,15,16), 2, 2, Col));
    /// ```
    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix) {
        assert_eq!(self.row, self.col);
        let r = self.row;
        let l = (self.row / 2) as i32;

        let mut v1 = vec![0f64; (l*l) as usize];
        let mut v2 = vec![0f64; l as usize * (r - l as usize)];
        let mut v3 = vec![0f64; (r - l as usize) * l as usize];
        let mut v4 = vec![0f64; (r - l as usize) * (r - l as usize)];

        for k in 0 .. (r*r) {
            let (quot, rem) = quot_rem(k, r);
            let r1 = r as i32;
            match (quot, rem) {
                (i, j) if (i < l) && (j < l) => {
                    let idx = (i * l + j) as usize;
                    v1[idx] = self.data[k];
                },
                (i, j) if (i < l) && (j >= l) => {
                    let idx = (i * (r1 - l) + j - l) as usize;
                    v2[idx] = self.data[k];
                },
                (i, j) if (i >= l) && (j < l) => {
                    let idx = ((i - l) * l + j) as usize;
                    v3[idx] = self.data[k];
                },
                (i, j) if (i >= l) && (j >= l) => {
                    let idx = ((i - l) * (r1 - l) + j - l) as usize;
                    v4[idx] = self.data[k];
                },
                _ => (),
            }
        }

        match self.shape {
            Row => {
                let m1 = matrix(v1, l as usize, l as usize, Row);
                let m2 = matrix(v2, l as usize, r - l as usize, Row);
                let m3 = matrix(v3, r - l as usize, l as usize, Row);
                let m4 = matrix(v4, r - l as usize, r - l as usize, Row);

                (m1, m2, m3, m4) // Row direction
            },
            Col => {
                let m1 = matrix(v1, l as usize, l as usize, Col);
                let m2 = matrix(v2, r - l as usize, l as usize, Col);
                let m3 = matrix(v3, l as usize, r - l as usize, Col);
                let m4 = matrix(v4, r - l as usize, r - l as usize, Col);

                (m1, m2, m3, m4) // Col direction
            }
        }
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
        let (p, l, u) = self.plu();
        if u.diag().into_iter().any(|x| x == 0f64) {
            None
        } else {
            Some(inv_u(u) % inv_l(l) % p)
        }
    }
}

// =============================================================================
// Back-end Utils
// =============================================================================

#[allow(unused_comparisons)]
pub fn tab(s: &str) -> String {
    let l = s.len();
    let mut m: String = String::new();
    if (5 - l) >= 0 {
        for _i in 0 .. (5 - l) {
            m += " ";
        }
    } else {
        m += " ";
    }
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
    let r1 = m1.row;
    let c2 = m4.col;

    let v1 = m1.data;
    let v2 = m2.data;
    let v3 = m3.data;
    let v4 = m4.data;

    let r = r1 + c2;
    let l = (r / 2) as i32;

    let mut v = vec![0f64; r*r];

    for k in 0 .. (r*r) {
        let (quot, rem) = quot_rem(k, r);
        let r1 = r as i32;
        match (quot, rem) {
            (i, j) if (i < l) && (j < l) => {
                let idx = (i * l + j) as usize;
                v[k] = v1[idx];
            },
            (i, j) if (i < l) && (j >= l) => {
                let idx = (i * (r1 - l) + j - l) as usize;
                v[k] = v2[idx];
            },
            (i, j) if (i >= l) && (j < l) => {
                let idx = ((i - l) * l + j) as usize;
                v[k] = v3[idx];
            },
            (i, j) if (i >= l) && (j >= l) => {
                let idx = ((i - l) * (r1 - l) + j - l) as usize;
                v[k] = v4[idx];
            },
            _ => (),
        }
    }

    matrix(v, r, r, m1.shape)
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
    let n = l.clone().data;
    let mut m = l.clone().data;

    match l.row {
        1 => l,
        2 => {
            m[2] = -n[2];
            matrix(m, 2, 2, Row)
        },
        _ => {
            let (l1, l2, l3, l4) = l.block();

            let m1 = inv_l(l1);
            let m2 = l2;
            let m4 = inv_l(l4);
            let m3 = ((m4.clone() % l3) % m1.clone()) * (-1f64);

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

    if w.shape == Col {
        w = w.change_shape();
    }

    let mut m = w.clone().data;
    let n = w.clone().data;

    match u.row {
        1 => {
            m[0] = 1f64 / n[0];
            matrix(m, 1, 1, Row)
        },
        2 => {
            let a = n[0];
            let b = n[1];
            let c = n[3];
            let d = a * c;

            m[0] = 1f64 / a;
            m[1] = -b / d;
            m[3] = 1f64 / c;

            matrix(m, 2, 2, Row)
        },
        _ => {
            let (u1, u2, u3, u4) = u.block();
            let m1 = inv_u(u1);
            let m3 = u3;
            let m4 = inv_u(u4);
            let m2 = (m1.clone() % u2 % m4.clone()) * (-1f64);

            combine(m1, m2, m3, m4)
        }
    }
}

/// Generate Pivot Matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = c!(0,1);
/// assert_eq!(pivot(a), matrix(c!(0,1,1,0),2,2,Row));
/// ```
pub fn pivot(v: Vector) -> Matrix {
    let c = v.len();
    let rank1 = v.rank();

    let mut data = vec![0f64; c*c];

    let mut j: usize = 0;

    for i in rank1.into_iter() {
        let idx = i * c + j;
        data[idx] = 1f64;
        j += 1;
    }
    matrix(data, c, c, Row)
}

/// Pivot to permutation
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let p = pivot(c!(1,3,2)); // [0 1 0; 0 0 1; 1 0 0]
/// let pm = to_perm(p);
/// assert_eq!(pm, c!(1,3,2).rank().into_iter().map(|x:usize| x as f64).collect::<Vec<f64>>());
/// ```
pub fn to_perm(p: Matrix) -> Vec<f64> {
    let mut m = p.clone();
    if m.shape != Col {
        m = m.change_shape();
    }
    let v = m.data;
    let mut container: Vec<f64> = Vec::new();
    for (i, x) in v.into_iter().enumerate() {
        if x == 1f64 {
            container.push((i % p.col) as f64);
        }
    }
    return container;
}

pub fn det_p(p: Matrix) -> Matrix {
    unimplemented!()
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