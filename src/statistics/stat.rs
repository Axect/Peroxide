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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
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
//!
//! ## Confusion Matrix
//!
//! * `ConfusionMatrix` is a struct to calculate confusion matrix
//! * The reference is [here](https://en.wikipedia.org/wiki/Confusion_matrix)
//!
//! ### Example
//!
//! ```rust
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let y     = vec![1usize, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0];
//!     let y_hat = vec![0usize, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0];
//!     let true_val = 1usize;
//!
//!     let cm = ConfusionMatrix::new(&y, &y_hat, true_val);
//!     cm.print();
//!     //         c[0]    c[1]
//!     // r[0]       6       2
//!     // r[1]       1       3
//!
//!     // to matrix
//!     let cm_mat = cm.to_matrix();
//!     
//!     // Calculate accuracy
//!     cm.ACC().print(); // 0.75
//!
//!     // Calculate TPR (Sensitivity or Recall)
//!     cm.TPR().print(); // 0.6666....
//!
//!     // Calculate some metrics
//!     let metrics = cm.calc_metrics(&[ACC, TPR, TNR, F1]);
//!
//!     // Print some metrics
//!     cm.summary(&[ACC, TPR, TNR, F1]);
//! }
//! ```

use std::fmt;

use self::QType::*;
//use crate::structure::dataframe::*;
use crate::structure::matrix::*;
use crate::traits::fp::FPVector;
use order_stat::kth_by;

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

impl Statistics for Vec<f64> {
    type Array = Vec<f64>;
    type Value = f64;

    /// Mean
    ///
    /// Uses welfords online algorithm for numerically stable computation.
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     assert_eq!(a.mean(), 3.0);
    /// }
    /// ```
    fn mean(&self) -> f64 {
        let mut xn = 0f64;
        let mut n = 0f64;

        for x in self.iter() {
            n += 1f64;
            xn += (x - xn) / n;
        }
        xn
    }

    /// Variance
    ///
    /// Uses welfords online algorithm for numerically stable computation.
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     assert_eq!(a.var(), 2.5);
    /// }
    /// ```
    fn var(&self) -> f64 {
        let mut xn = 0f64;
        let mut n = 0f64;
        let mut m2n: f64 = 0f64;

        for x in self.iter() {
            n += 1f64;
            let diff_1 = x - xn;
            xn += diff_1 / n;
            m2n += diff_1 * (x - xn);
        }
        assert_ne!(n, 1f64);
        m2n / (n - 1f64)
    }

    /// Standard Deviation
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3);
    ///     assert!(nearly_eq(a.sd(), 1f64)); // Floating Number Error
    /// }
    /// ```
    fn sd(&self) -> f64 {
        self.var().sqrt()
    }

    fn cov(&self) -> Vec<f64> {
        unimplemented!()
    }
    fn cor(&self) -> Vec<f64> {
        unimplemented!()
    }
}

impl Statistics for Vec<f32> {
    type Array = Vec<f32>;
    type Value = f32;

    /// Mean
    ///
    /// Uses welfords online algorithm for numerically stable computation.
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     assert_eq!(a.mean(), 3.0);
    /// }
    /// ```
    fn mean(&self) -> f32 {
        let mut xn = 0f32;
        let mut n = 0f32;

        for x in self.iter() {
            n += 1f32;
            xn += (x - xn) / n;
        }
        xn
    }

    /// Variance
    ///
    /// Uses welfords online algorithm for numerically stable computation.
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     assert_eq!(a.var(), 2.5);
    /// }
    /// ```
    fn var(&self) -> f32 {
        let mut xn = 0f32;
        let mut n = 0f32;
        let mut m2n: f32 = 0f32;

        for x in self.iter() {
            n += 1f32;
            let diff_1 = x - xn;
            xn += diff_1 / n;
            m2n += diff_1 * (x - xn);
        }
        assert_ne!(n, 1f32);
        m2n / (n - 1f32)
    }

    /// Standard Deviation
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3);
    ///     assert!(nearly_eq(a.sd(), 1f64)); // Floating Number Error
    /// }
    /// ```
    fn sd(&self) -> f32 {
        self.var().sqrt()
    }

    fn cov(&self) -> Vec<f32> {
        unimplemented!()
    }
    fn cor(&self) -> Vec<f32> {
        unimplemented!()
    }
}

impl Statistics for Matrix {
    type Array = Matrix;
    type Value = Vec<f64>;

    /// Column Mean
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let m = matrix(c!(1,3,3,1), 2, 2, Col);
    ///     assert_eq!(m.mean(), c!(2,2));
    /// }
    /// ```
    fn mean(&self) -> Vec<f64> {
        let mut container: Vec<f64> = Vec::new();
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
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let m = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    ///     assert!(nearly_eq(m.var()[0], 1));
    /// }
    /// ```
    fn var(&self) -> Vec<f64> {
        let mut container: Vec<f64> = Vec::new();
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
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let m = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    ///     assert!(nearly_eq(m.sd()[0], 1));
    /// }
    /// ```
    fn sd(&self) -> Vec<f64> {
        let mut container: Vec<f64> = Vec::new();
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
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let m = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    ///     println!("{}", m.cov());
    ///
    ///     //         c[0]    c[1]
    ///     // r[0]  1.0000 -1.0000
    ///     // r[1] -1.0000  1.0000
    /// }
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

//impl Statistics for DataFrame {
//    type Array = Matrix;
//    type Value = Self;
//
//    fn mean(&self) -> Self::Value {
//        let mut df = DataFrame::with_header(self.headers().map(|x| x.as_str()).collect());
//        for k in self.headers() {
//            df[k] = vec![self[k].mean()];
//        }
//        df
//    }
//
//    fn var(&self) -> Self::Value {
//        let mut df = DataFrame::with_header(self.headers().map(|x| x.as_str()).collect());
//        for k in self.headers() {
//            df[k] = vec![self[k].var()];
//        }
//        df
//    }
//
//    fn sd(&self) -> Self::Value {
//        let mut df = DataFrame::with_header(self.headers().map(|x| x.as_str()).collect());
//        for k in self.headers() {
//            df[k] = vec![self[k].sd()];
//        }
//        df
//    }
//
//    fn cov(&self) -> Self::Array {
//        let m: Matrix = self.into();
//        m.cov()
//    }
//
//    fn cor(&self) -> Self::Array {
//        self.to_matrix().cor()
//    }
//}

/// Covariance (to Value)
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let v1 = c!(1,2,3);
///     let v2 = c!(3,2,1);
///     assert!(nearly_eq(cov(&v1, &v2), -1f64));
/// }
/// ```
pub fn cov(v1: &Vec<f64>, v2: &Vec<f64>) -> f64 {
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
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = c!(1,2,3);
///     let b = c!(3,2,1);
///     assert!(nearly_eq(cor(&a, &b),-1));
/// }
/// ```
pub fn cor(v1: &Vec<f64>, v2: &Vec<f64>) -> f64 {
    cov(v1, v2) / (v1.sd() * v2.sd())
}

/// R like linear regression
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a: Matrix = c!(1,2,3,4,5).into();
///     let b: Matrix = &a + &Normal(0,1).sample(5).into();
///     lm(&a, &b).print();
///
///     //        c[0]
///     // r[0] 0.7219
///     // r[1] 0.8058
/// }
/// ```
pub fn lm(input: &Matrix, target: &Matrix) -> Matrix {
    let mut ones = vec![1f64; input.row * input.col];
    ones.extend(&input.data);
    let x = matrix(ones, input.row, input.col + 1, input.shape);
    &x.pseudo_inv() * target
}

// =============================================================================
// Ordered Statistics (Use `order-stat`)
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

impl OrderedStat for Vec<f64> {
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
        for i in 0..q.len() {
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
            } else if q == 1f64 {
                l - 1
            } else if q - (k as f64) * p > 0f64 {
                k
            } else {
                k - 1
            };
            *kth_by(v, k, |x, y| x.partial_cmp(y).unwrap())
        }
        Type2 => {
            if q == 0f64 {
                let k = 0;
                *kth_by(v, k, |x, y| x.partial_cmp(y).unwrap())
            } else if q == 1f64 {
                let k = l - 1;
                *kth_by(v, k, |x, y| x.partial_cmp(y).unwrap())
            } else if q - (k as f64) * p > 0f64 {
                *kth_by(v, k, |x, y| x.partial_cmp(y).unwrap())
            } else {
                let prev = *kth_by(v, k - 1, |x, y| x.partial_cmp(y).unwrap());
                let next = *kth_by(v, k, |x, y| x.partial_cmp(y).unwrap());
                (prev + next) / 2f64
            }
        }
        _ => unimplemented!(),
    }
}

pub fn quantile(v: &Vec<f64>, qtype: QType) -> Vec<f64> {
    let q_vec = vec![0.0, 0.25, 0.5, 0.75, 1.0];
    v.quantiles(q_vec, qtype)
}

// =============================================================================
// Confusion Matrix
// =============================================================================
/// Confusion Matrix
/// 
/// * `TP` : True Positive
/// * `TN` : True Negative
/// * `FP` : False Positive
/// * `FN` : False Negative
///
/// # Examples
/// ```
/// use peroxide::fuga::*;
///
/// fn main() {
///     let y           = vec![1usize, 1, 1, 0, 0, 0];
///     let y_hat       = vec![1usize, 0, 1, 0, 0, 1];
///     let true_val    = 1usize;
///
///     // Create Confusion Matrix
///     let cm = ConfusionMatrix::new(&y, &y_hat, true_val);
///
///     // Print
///     cm.print();
///     //        c[0]  c[1]
///     // r[0]  2.0000  1.0000
///     // r[1]  1.0000  2.0000
///
///     // To Matrix
///     let cm_mat = cm.to_matrix();
///
///     // Calculate Accuracy
///     let acc = cm.ACC();
///
///     // Calculate for some metrics
///     let metrics = cm.calc_metrics(&[ACC, TPR, FPR, F1]);
///
///     // Print summary for some metrics
///     cm.summary(&[ACC, TPR, FPR, F1]);
/// }
/// ```
#[allow(non_snake_case)]
#[derive(Debug, Clone, PartialEq)]
pub struct ConfusionMatrix {
    pub TP: usize,
    pub TN: usize,
    pub FP: usize,
    pub FN: usize,
}

impl ConfusionMatrix {
    /// Create Confusion Matrix
    /// 
    /// # Examples
    /// ```
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let y           = vec![1usize, 1, 1, 0, 0, 0];
    ///     let y_hat       = vec![1usize, 0, 1, 0, 0, 1];
    ///     let true_val    = 1usize;
    ///
    ///     let true_val = 1usize;
    ///     let cm = ConfusionMatrix::new(&y, &y_hat, true_val);
    ///     cm.print();
    ///     //        c[0]  c[1]
    ///     // r[0]  2.0000  1.0000
    ///     // r[1]  1.0000  2.0000
    /// }
    /// ```
    #[allow(non_snake_case)]
    pub fn new<T: PartialEq + Clone + Copy>(y: &Vec<T>, y_hat: &Vec<T>, true_val: T) -> Self {
        let mut TP = 0;
        let mut TN = 0;
        let mut FP = 0;
        let mut FN = 0;

        for (&y, &y_hat) in y.iter().zip(y_hat.iter()) {
            if y == true_val && y_hat == true_val {
                TP += 1;
            } else if y != true_val && y_hat != true_val {
                TN += 1;
            } else if y != true_val && y_hat == true_val {
                FP += 1;
            } else if y == true_val && y_hat != true_val {
                FN += 1;
            }
        }

        Self {
            TP,
            TN,
            FP,
            FN,
        }
    }

    /// Condition Positive
    #[allow(non_snake_case)]
    pub fn P(&self) -> usize {
        self.TP + self.FN
    }

    /// Condition Negative
    #[allow(non_snake_case)]
    pub fn N(&self) -> usize {
        self.TN + self.FP
    }

    /// True Positive Rate (Sensitivity, Recall, Hit-rate)
    #[allow(non_snake_case)]
    pub fn TPR(&self) -> f64 {
        self.TP as f64 / (self.TP + self.FN) as f64
    }

    /// True Negative Rate (Specificity, Selectivity)
    #[allow(non_snake_case)]
    pub fn TNR(&self) -> f64 {
        self.TN as f64 / (self.TN + self.FP) as f64
    }

    /// Positive Predictive Value (Precision)
    #[allow(non_snake_case)]
    pub fn PPV(&self) -> f64 {
        self.TP as f64 / (self.TP + self.FP) as f64
    }

    /// Negative Predictive Value
    #[allow(non_snake_case)]
    pub fn NPV(&self) -> f64 {
        self.TN as f64 / (self.TN + self.FN) as f64
    }

    /// False Negative Rate (Miss-rate)
    #[allow(non_snake_case)]
    pub fn FNR(&self) -> f64 {
        self.FN as f64 / (self.FN + self.TP) as f64
    }

    /// False Positive Rate (Fall-out)
    #[allow(non_snake_case)]
    pub fn FPR(&self) -> f64 {
        self.FP as f64 / (self.FP + self.TN) as f64
    }

    /// False Discovery Rate
    #[allow(non_snake_case)]
    pub fn FDR(&self) -> f64 {
        self.FP as f64 / (self.FP + self.TP) as f64
    }

    /// False Omission Rate
    #[allow(non_snake_case)]
    pub fn FOR(&self) -> f64 {
        self.FN as f64 / (self.FN + self.TN) as f64
    }

    /// Positive Likelihood Ratio
    #[allow(non_snake_case)]
    pub fn LR_plus(&self) -> f64 {
        self.TPR() / self.FPR()
    }

    /// Negative Likelihood Ratio
    #[allow(non_snake_case)]
    pub fn LR_minus(&self) -> f64 {
        self.FNR() / self.TNR()
    }

    /// Prevalence Threshold
    #[allow(non_snake_case)]
    pub fn PT(&self) -> f64 {
        self.FPR().sqrt() / (self.TPR().sqrt() + self.FPR().sqrt())
    }

    /// Threat Score (Critical Success Index)
    #[allow(non_snake_case)]
    pub fn TS(&self) -> f64 {
        self.TP as f64 / (self.TP + self.FP + self.FN) as f64
    }

    /// Prevalence
    #[allow(non_snake_case)]
    pub fn prevalence(&self) -> f64 {
        self.P() as f64 / (self.P() + self.N()) as f64
    }

    /// Accuracy
    #[allow(non_snake_case)]
    pub fn ACC(&self) -> f64 {
        (self.TP + self.TN) as f64 / (self.P() + self.N()) as f64
    }

    /// Balanced Accuracy
    #[allow(non_snake_case)]
    pub fn BA(&self) -> f64 {
        (self.TPR() + self.TNR()) / 2f64
    }

    /// F1 Score
    #[allow(non_snake_case)]
    pub fn F1(&self) -> f64 {
        2f64 * self.TP as f64 / (2f64 * self.TP as f64 + self.FP as f64 + self.FN as f64)
    }

    /// Matthews Correlation Coefficient (Phi Coefficient)
    #[allow(non_snake_case)]
    pub fn MCC(&self) -> f64 {
        let a = self.TP as f64;
        let b = self.FP as f64;
        let c = self.FN as f64;
        let d = self.TN as f64;
        (a * d - b * c) / ((a + b) * (a + c) * (d + b) * (d + c)).sqrt()
    }

    /// Fowlkes-Mallows Index
    #[allow(non_snake_case)]
    pub fn FM(&self) -> f64 {
        (self.PPV() * self.TPR()).sqrt()
    }

    /// Bookmaker Informedness (Informedness)
    #[allow(non_snake_case)]
    pub fn BM(&self) -> f64 {
        self.TPR() + self.TNR() - 1f64
    }

    /// Markedness (deltaP)
    #[allow(non_snake_case)]
    pub fn MK(&self) -> f64 {
        self.PPV() + self.NPV() - 1f64
    }

    /// Diagnostic Odds Ratio
    #[allow(non_snake_case)]
    pub fn DOR(&self) -> f64 {
        self.LR_plus() / self.LR_minus()
    }

    /// To Matrix
    pub fn to_matrix(&self) -> Matrix {
        let mut m = matrix(vec![0f64; 4], 2, 2, Row);
        m[(0,0)] = self.TP as f64;
        m[(0,1)] = self.FP as f64;
        m[(1,0)] = self.FN as f64;
        m[(1,1)] = self.TN as f64;
        m
    }

    /// Calculate a specific metric
    pub fn calc_metric(&self, metric: Metric) -> f64 {
        match metric {
            Metric::TPR => self.TPR(),
            Metric::TNR => self.TNR(),
            Metric::PPV => self.PPV(),
            Metric::NPV => self.NPV(),
            Metric::FNR => self.FNR(),
            Metric::FPR => self.FPR(),
            Metric::FDR => self.FDR(),
            Metric::FOR => self.FOR(),
            Metric::LR_plus => self.LR_plus(),
            Metric::LR_minus => self.LR_minus(),
            Metric::PT => self.PT(),
            Metric::TS => self.TS(),
            Metric::prevalence => self.prevalence(),
            Metric::ACC => self.ACC(),
            Metric::BA => self.BA(),
            Metric::F1 => self.F1(),
            Metric::MCC => self.MCC(),
            Metric::FM => self.FM(),
            Metric::BM => self.BM(),
            Metric::MK => self.MK(),
            Metric::DOR => self.DOR(),
        }
    }

    /// Calculate for some metrics
    pub fn calc_metrics(&self, metrics: &[Metric]) -> Vec<f64> {
        metrics.iter().map(|m| self.calc_metric(*m)).collect()
    }

    /// Summarize some metrics
    pub fn summary(&self, metrics: &[Metric]) {
        let width = metrics.iter().fold(0, |acc, m| {
            if m.to_string().len() > acc {
                m.to_string().len()
            } else {
                acc
            }
        });
        println!("============================================================");
        println!("Summary of metrics");
        println!("============================================================");
        for m in metrics {
            println!("{:width$}:\t{:.4}", m.to_string(), self.calc_metric(*m));
        }
        println!("============================================================");
    }
}

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy)]
pub enum Metric {
    TPR,
    TNR,
    PPV,
    NPV,
    FNR,
    FPR,
    FDR,
    FOR,
    LR_plus,
    LR_minus,
    PT,
    TS,
    prevalence,
    ACC,
    BA,
    F1,
    MCC,
    FM,
    BM,
    MK,
    DOR,
}

impl fmt::Display for Metric {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Metric::TPR => write!(f, "TPR"),
            Metric::TNR => write!(f, "TNR"),
            Metric::PPV => write!(f, "PPV"),
            Metric::NPV => write!(f, "NPV"),
            Metric::FNR => write!(f, "FNR"),
            Metric::FPR => write!(f, "FPR"),
            Metric::FDR => write!(f, "FDR"),
            Metric::FOR => write!(f, "FOR"),
            Metric::LR_plus => write!(f, "LR+"),
            Metric::LR_minus => write!(f, "LR-"),
            Metric::PT => write!(f, "PT"),
            Metric::TS => write!(f, "TS"),
            Metric::prevalence => write!(f, "prevalence"),
            Metric::ACC => write!(f, "ACC"),
            Metric::BA => write!(f, "BA"),
            Metric::F1 => write!(f, "F1"),
            Metric::MCC => write!(f, "MCC"),
            Metric::FM => write!(f, "FM"),
            Metric::BM => write!(f, "BM"),
            Metric::MK => write!(f, "MK"),
            Metric::DOR => write!(f, "DOR"),
        }
    }
}
