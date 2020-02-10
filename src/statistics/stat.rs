//! Basic statistics
//!
//! ## `Statistics` trait
//!
//! * To make generic code, there is `Statistics` trait
//!     * `mean`: just mean
//!     * `var` : variance
//!     * `sd` : standard deviation (R-like notation)
//!     * `cov` : covariance
//!     * `cor` : correlation coefficient
//!     ```rust
//!     pub trait Statistics {
//!         type Array;
//!         type Value;
//!
//!         fn mean(&self) -> Self::Value;
//!         fn var(&self) -> Self::Value;
//!         fn sd(&self) -> Self::Value;
//!         fn cov(&self) -> Self::Array;
//!         fn cor(&self) -> Self::Array;
//!     }
//!     ```
//!
//! ### For `Vec<f64>`
//!
//! * Caution: For `Vec<f64>`, `cov` & `cor` are unimplemented (those for `Matrix`)
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4,5);
//!         a.mean().print(); // 3
//!         a.var().print();  // 2.5
//!         a.sd().print();   // 1.5811388300841898
//!     }
//!     ```
//!
//! * But there are other functions to calculate `cov` & `cor`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let v1 = c!(1,2,3);
//!         let v2 = c!(3,2,1);
//!
//!         cov(&v1, &v2).print(); // -0.9999999999999998
//!         cor(&v1, &v2).print(); // -0.9999999999999993
//!     }
//!     ```
//!
//! ### For `Matrix`
//!
//! * For `Matrix`, `mean, var, sd` means column operations
//! * `cov` means covariance matrix & `cor` means also correlation coefficient matrix
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let m = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
//!
//!         m.mean().print(); // [2, 2]
//!         m.var().print();  // [1.0000, 1.0000]
//!         m.sd().print();   // [1.0000, 1.0000]
//!
//!         m.cov().print();
//!         //         c[0]    c[1]
//!         // r[0]  1.0000 -1.0000
//!         // r[1] -1.0000  1.0000
//!
//!         m.cor().print();
//!         //         c[0]    c[1]
//!         // r[0]       1 -1.0000
//!         // r[1] -1.0000       1
//!     }
//!     ```
//! 
//! ### For `DataFrame`
//! 
//! * Similar to Matrix but, `Value` is `DataFrame`
//! * `cov` means covariance matrix.
//! 
//! ```rust
//! extern crate peroxide;
//! use peroxide::*;
//! 
//! fn main() {
//!     #[cfg(feature = "dataframe")]
//!     {
//!         let mut m = DataFrame::with_header(vec!["x", "y"]);
//!         m["x"] = c!(1,2,3);
//!         m["y"] = c!(3,2,1);
//!         
//!         m.cov().print();
//!         //         c[0]    c[1]
//!         // r[0]  1.0000 -1.0000
//!         // r[1] -1.0000  1.0000
//!     }
//! }
//! ```

use structure::matrix::*;
use structure::vector::*;
#[cfg(feature = "dataframe")]
use structure::dataframe::*;
use order_stat::{kth_by};
use self::QType::*;

/// Statistics Trait
///
/// It contains `mean`, `var`, `sd`, `cov`
pub trait Statistics {
    type Array;
    type Value;

    fn mean(&self) -> Self::Value;
    fn var(&self) -> Self::Value;
    fn sd(&self) -> Self::Value;
    fn cov(&self) -> Self::Array;
    fn cor(&self) -> Self::Array;
}

impl Statistics for Vector {
    type Array = Vector;
    type Value = f64;

    /// Mean
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// assert_eq!(a.mean(), 3.0);
    /// ```
    fn mean(&self) -> f64 {
        self.reduce(0f64, |x, y| x + y) / (self.len() as f64)
    }

    /// Variance
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// assert_eq!(a.var(), 2.5);
    /// ```
    fn var(&self) -> f64 {
        let mut ss = 0f64;
        let mut s = 0f64;
        let mut l = 0f64;

        for x in self.into_iter() {
            ss += x.powf(2f64);
            s += *x;
            l += 1f64;
        }
        assert_ne!(l, 1f64);
        (ss / l - (s / l).powf(2f64)) * l / (l - 1f64)
    }

    /// Standard Deviation
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3);
    /// assert!(nearly_eq(a.sd(), 1f64)); // Floating Number Error
    /// ```
    fn sd(&self) -> f64 {
        self.var().sqrt()
    }

    fn cov(&self) -> Vector {
        unimplemented!()
    }
    fn cor(&self) -> Vector {
        unimplemented!()
    }
}

impl Statistics for Matrix {
    type Array = Matrix;
    type Value = Vector;

    /// Column Mean
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let m = matrix(c!(1,3,3,1), 2, 2, Col);
    /// assert_eq!(m.mean(), c!(2,2));
    /// ```
    fn mean(&self) -> Vector {
        let mut container: Vector = Vec::new();
        let c = self.col;

        for i in 0..c {
            container.push(self.col(i).mean());
        }
        container
    }

    /// Column variance
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let m = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// assert!(nearly_eq(m.var()[0], 1));
    /// ```
    fn var(&self) -> Vector {
        let mut container: Vector = Vec::new();
        let c = self.col;

        for i in 0..c {
            container.push(self.col(i).var());
        }
        container
    }

    /// Column Standard Deviation
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let m = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// assert!(nearly_eq(m.sd()[0], 1));
    /// ```
    fn sd(&self) -> Vector {
        let mut container: Vector = Vec::new();
        let c = self.col;

        for i in 0..c {
            container.push(self.col(i).sd());
        }
        container
    }

    /// Covariance Matrix (Column based)
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let m = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// println!("{}", m.cov());
    ///
    /// //         c[0]    c[1]
    /// // r[0]  1.0000 -1.0000
    /// // r[1] -1.0000  1.0000
    /// ```
    fn cov(&self) -> Self {
        let c = self.col;

        let mut m: Self = matrix(vec![0f64; c * c], c, c, self.shape);

        for i in 0..c {
            let m1 = self.col(i);
            for j in 0..c {
                let m2 = self.col(j);
                m[(i, j)] = cov(&m1, &m2);
            }
        }
        m
    }

    fn cor(&self) -> Self {
        let c = self.col;

        let mut m: Self = matrix(vec![0f64; c * c], c, c, self.shape);

        for i in 0..c {
            let m1 = self.col(i);
            for j in 0..c {
                let m2 = self.col(j);
                m[(i, j)] = cor(&m1, &m2);
            }
        }
        m
    }
}

#[cfg(feature = "dataframe")]
impl Statistics for DataFrame {
    type Array = Matrix;
    type Value = Self;

    fn mean(&self) -> Self::Value {
        let mut df = DataFrame::with_header(self.headers().map(|x| x.as_str()).collect());
        for k in self.headers() {
            df[k] = vec![self[k].mean()];
        }
        df
    }

    fn var(&self) -> Self::Value {
        let mut df = DataFrame::with_header(self.headers().map(|x| x.as_str()).collect());
        for k in self.headers() {
            df[k] = vec![self[k].var()];
        }
        df
    }

    fn sd(&self) -> Self::Value {
        let mut df = DataFrame::with_header(self.headers().map(|x| x.as_str()).collect());
        for k in self.headers() {
            df[k] = vec![self[k].sd()];
        }
        df
    }

    fn cov(&self) -> Self::Array {
        self.to_matrix().cov()
    }

    fn cor(&self) -> Self::Array {
        self.to_matrix().cor()
    }
}

/// Covariance (to Value)
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let v1 = c!(1,2,3);
/// let v2 = c!(3,2,1);
/// assert!(nearly_eq(cov(&v1, &v2), -1f64));
/// ```
pub fn cov(v1: &Vector, v2: &Vector) -> f64 {
    let mut ss = 0f64;
    let mut sx = 0f64;
    let mut sy = 0f64;
    let mut l = 0f64;

    for (x, y) in v1.into_iter().zip(v2) {
        ss += x * y;
        sx += *x;
        sy += *y;
        l += 1f64;
    }
    assert_ne!(l, 1f64);
    (ss / l - (sx * sy) / (l.powf(2f64))) * l / (l - 1f64)
}

/// Pearson's correlation coefficient
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = c!(1,2,3);
/// let b = c!(3,2,1);
/// assert!(nearly_eq(cor(&a, &b),-1));
/// ```
pub fn cor(v1: &Vector, v2: &Vector) -> f64 {
    cov(v1, v2) / (v1.sd() * v2.sd())
}

/// R like linear regression
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = c!(1,2,3,4,5).to_matrix();
/// let b = &a + &Normal(0,1).sample(5).to_matrix();
/// lm(&a, &b).print();
///
/// //        c[0]
/// // r[0] 0.7219
/// // r[1] 0.8058
/// ```
pub fn lm(input: &Matrix, target: &Matrix) -> Matrix {
    let mut ones = vec![1f64; input.row * input.col];
    ones.extend(&input.data);
    let x = matrix(ones, input.row, input.col + 1, input.shape);
    &x.pseudo_inv().unwrap() * target
}

// =============================================================================
// Ordered Statistics (Use `order-stat`
// =============================================================================
/// Trait for Ordered Statistics
///
/// * `median`
/// * `quantile`
pub trait OrderedStat {
    type Array;
    type Value;

    fn median(&self) -> Self::Value;
    fn quantile(&self, q: f64, qtype: QType) -> Self::Value;
    fn quantiles(&self, q: Vec<f64>, qtype: QType) -> Self::Array;
}

/// R Quantile Type enums
#[derive(Debug, Copy, Clone)]
pub enum QType {
    Type1,
    Type2,
    Type3,
    Type4,
    Type5,
    Type6,
    Type7,
    Type8,
    Type9,
}

impl OrderedStat for Vector {
    type Array = Self;
    type Value = f64;

    fn median(&self) -> Self::Value {
        self.quantile(0.5, Type2)
    }

    fn quantile(&self, q: f64, qtype: QType) -> Self::Value {
        let mut m = self.clone();
        quantile_mut(&mut m, q, qtype)
    }

    fn quantiles(&self, q: Vec<f64>, qtype: QType) -> Self::Array {
        let mut v = vec![0f64; q.len()];
        let mut m = self.clone();
        for i in 0 .. q.len() {
            v[i] = quantile_mut(&mut m, q[i], qtype);
        }
        v
    }
}

fn quantile_mut(v: &mut [f64], q: f64, t: QType) -> f64 {
    let l = v.len();
    let p = 1f64 / (l as f64);
    let k = (q / p) as usize;
    match t {
        Type1 => {
            let k = if q == 0f64 {
                0
            } else if q - (k as f64) * p > 0f64 {
                k
            } else {
                k - 1
            };
            *kth_by(v, k, |x,y| x.partial_cmp(y).unwrap())
        }
        Type2 => {
            if q - (k as f64) * p > 0f64 {
                *kth_by(v, k, |x,y| x.partial_cmp(y).unwrap())
            } else if q == 0f64 {
                let k = 0;
                *kth_by(v, k, |x,y| x.partial_cmp(y).unwrap())
            } else if q == 1f64 {
                let k = l - 1;
                *kth_by(v, k, |x,y| x.partial_cmp(y).unwrap())
            } else {
                let prev = *kth_by(v, k-1, |x,y| x.partial_cmp(y).unwrap());
                let next = *kth_by(v, k, |x,y| x.partial_cmp(y).unwrap());
                (prev + next) / 2f64
            }
        }
        _ => unimplemented!()
    }
}

pub fn quantile(v: &Vec<f64>, qtype: QType) -> Vec<f64> {
    let q_vec = vec![0.0, 0.25, 0.5, 0.75, 1.0];
    v.quantiles(q_vec, qtype)
}
