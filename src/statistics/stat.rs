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

use structure::matrix::*;
use structure::vector::*;

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
        sx += x;
        sy += y;
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
/// let b = a.clone() + Normal(0,1).sample(5).to_matrix();
/// lm!(b ~ a).print();
///
/// //        c[0]
/// // r[0] 0.7219
/// // r[1] 0.8058
/// ```
pub fn lm(input: &Matrix, target: &Matrix) -> Matrix {
    let x_temp = input.clone();
    let mut ones = vec![1f64; x_temp.row * x_temp.col];
    ones.extend(&x_temp.data);
    let x = matrix(ones, x_temp.row, x_temp.col + 1, x_temp.shape);
    let t = target.clone();
    x.pseudo_inv().unwrap() * t
}
