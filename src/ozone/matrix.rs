use std::convert;
use std::fmt;
use std::iter::FromIterator;
use std::ops::{Add, Neg, Sub, Mul};
pub use self::Shape::{Row, Col};

/// To select matrices' binding.
/// 
/// Row - Row binding
/// Col - Column binding
///
/// # Examples
/// ```
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
        let to_disp = match self {
            Row => "Row",
            Col => "Col",
        };
        write!(f, "{}", to_disp)
    }
}

/// R-like matrix structure
///
/// # Examples
///
/// ```
/// let a = Matrix {
///     data: vec![1f64,2,3,4],
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
pub trait Generic<T: convert::Into<f64>> {
    fn new(v: Vec<T>, x:usize, y:usize, shape: Shape) -> Matrix;
}

impl<T> Generic<T> for Matrix where T: convert::Into<f64> {
    /// Matrix generic constructor
    ///
    /// You can use any numeric type vector
    /// e.g. `u32`, `i32`, `i64`, ...
    ///
    /// # Examples
    /// ```
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

/// Pretty Print
impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

#[allow(dead_code)]
impl Matrix {
    /// Change Bindings
    ///
    /// `Row` -> `Col` or `Col` -> `Row`
    ///
    /// # Examples
    /// ```
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
    /// let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
    /// println!("{}", a) /// [[1,3],[2,4]]
    /// ```
    pub fn transpose(&self) -> Matrix {
        match self.shape {
            Row => Matrix::new(self.data.clone(), self.col, self.row, Col),
            Col => Matrix::new(self.data.clone(), self.col, self.row, Row)
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
impl Add for Matrix {
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

/// Negation of Matrix
///
/// # Examples
/// ```
/// let a = Matrix::new(vec![1,2,3,4],2,2,Row);
/// println!("{}", -a); // [[-1,-2],[-3,-4]]
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

impl Sub for Matrix {
    type Output = Matrix;

    fn sub(self, other: Matrix) -> Matrix {
        self.add(other.neg())
    }
}

impl Mul for Matrix {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
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
            self.change_shape().mul(other.change_shape())
        } else if self.shape == Col {
            self.change_shape().mul(other)
        } else {
            self.mul(other.change_shape())
        }
    }
}

//TODO: Implement iterator for Matrix
// pub struct MatrixIter {
//     matrix: Matrix,
//     index: usize,
// }

// impl IntoIterator for Matrix {
//     type Item = Matrix;
//     type IntoIter = MatrixIter;

//     fn into_iter(self) -> Self::IntoIter {
//         MatrixIter {
//             matrix: self,
//             index: 0,
//         }
//     }
// }

// impl Iterator for MatrixIter {
//     type Item = Matrix;

//     fn next(&mut self) -> Option<Self::Item> {
//         let mut result: Vec<f64> = Vec::new();
//         match self.matrix.shape {
//             Row => {
//                 if self.index >= self.matrix.row {
//                     return None;
//                 }
//                 result = self.matrix.data.clone().into_iter().skip(self.index * self.matrix.col).take(self.matrix.col).collect::<Vec<f64>>();
//                 self.index += 1;
//                 Some(
//                     Matrix::new(
//                         result,
//                         1,
//                         self.matrix.col,
//                         Row,
//                     )
//                 )
//             },
//             Col => {
//                 if self.index >= self.matrix.col {
//                     return None;
//                 }
//                 result = self.matrix.data.clone().into_iter().skip(self.index * self.matrix.row).take(self.matrix.row).collect::<Vec<f64>>();
//                 self.index += 1;
//                 Some(
//                     Matrix::new(
//                         result,
//                         self.matrix.row,
//                         1,
//                         Col,
//                     )
//                 )
//             }
//         }
//     }
// }

// impl FromIterator<Matrix> for Matrix {
//     fn from_iter<I: IntoIterator<Item=Matrix>(iter: I) -> Self {
//         let 
//     }
// }

// =============================================================================
// Functional Programming Tools (Hand-written)
// =============================================================================
pub trait FP {
    fn fmap<F>(&self, f: F) -> Matrix where F: Fn(f64) -> f64;
    fn reduce<F, T>(&self, init: T, f: F) -> f64 
        where F: Fn(f64, f64) -> f64,
              T: convert::Into<f64>;
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
