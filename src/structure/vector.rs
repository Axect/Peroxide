//! Extra tools for `Vec<f64>`
//!
//! ## Print `Vec<f64>`
//!
//! * There are two ways to print vector
//!     * Original way: `print!("{:?}", a);`
//!     * Peroxide way: `a.print();` - **Round-off to fourth digit**
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = vec![2f64.sqrt()];
//!         a.print(); // [1.4142]
//!     }
//!     ```
//!
//! ## Syntactic sugar for `Vec<f64>`
//!
//! * There is useful macro for `Vec<f64>`
//! * For `R`, there is `c`
//!
//!     ```R
//!     # R
//!     a = c(1,2,3,4)
//!     ```
//!
//! * For `Peroxide`, there is `c!`
//!
//!     ```rust
//!     // Rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!     }
//!     ```
//!
//! ## From ranges to Vector
//!
//! * For `R`, there is `seq` to declare sequence.
//!
//!     ```R
//!     # R
//!     a = seq(1, 4, 1)
//!     print(a)
//!     # [1] 1 2 3 4
//!     ```
//!
//! * For `peroxide`, there is `seq` to declare sequence.
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = seq(1, 4, 1);
//!         a.print();
//!         // [1, 2, 3, 4]
//!     }
//!     ```
//!
//! ## Vector Operation
//!
//! * There are some vector-wise operations
//!     * `add(&self, other: Vec<f64>) -> Vec<f64>`
//!     * `sub(&self, other: Vec<f64>) -> Vec<f64>`
//!     * `mul(&self, other: Vec<f64>) -> Vec<f64>`
//!     * `div(&self, other: Vec<f64>) -> Vec<f64>`
//!     * `dot(&self, other: Vec<f64>) -> f64`
//!     * `norm(&self) -> f64`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         let b = c!(4,3,2,1);
//!
//!         a.add(&b).print();
//!         a.sub(&b).print();
//!         a.mul(&b).print();
//!         a.div(&b).print();
//!         a.dot(&b).print();
//!         a.norm().print();
//!
//!         // [5, 5, 5, 5]
//!         // [-3, -1, 1, 3]
//!         // [4, 6, 6, 4]
//!         // [0.25, 0.6667, 1.5, 4]
//!         // 20
//!         // 5.477225575051661 // sqrt(30)
//!     }
//!     ```
//!
//! ## Concatenation
//!
//! There are two concatenation operations.
//!
//! * `cat(T, Vec<T>) -> Vec<f64>`
//! * `concat(Vec<T>, Vec<T>) -> Vec<T>`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         cat(0f64, a.clone()).print();
//!         // [0, 1, 2, 3, 4]
//!
//!         let b = c!(5,6,7,8);
//!         concat(a, b).print();
//!         // [1, 2, 3, 4, 5, 6, 7, 8]
//!     }
//!     ```
//!
//! ## Conversion to Matrix
//!
//! There are two ways to convert vector to matrix.
//!
//! * `to_matrix(&self) -> Matrix` : Vector to column matrix
//! * `transpose(&self) -> Matrix` : Vector to row matrix
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         let m_col = matrix(c!(1,2,3,4), 4, 1, Col); // (4,1) Matrix
//!         assert_eq!(a.to_matrix(), m_col);
//!
//!         let m_row = matrix(c!(1,2,3,4), 1, 4, Row); // (1,4) Matrix
//!         assert_eq!(a.transpose(), m_row);
//!     }
//!     ```
//!
//! # Functional Programming {#functional}
//!
//! ## FP for Vector
//!
//! * There are some functional programming tools for `Vec<f64>`
//!
//! ### fmap
//!
//! * `fmap` is syntactic sugar for `map`
//! * But different to original `map` - Only `f64 -> f64` allowed.
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!
//!         // Original rust
//!         a.clone()
//!             .into_iter()
//!             .map(|x| x + 1f64)
//!             .collect::<Vec<f64>>()
//!             .print();
//!             // [2, 3, 4, 5]
//!
//!         // fmap in Peroxide
//!         a.fmap(|x| x + 1f64).print();
//!         // [2, 3, 4, 5]
//!     }
//!     ```
//!
//! ### reduce
//!
//! * `reduce` is syntactic sugar for `fold`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!
//!         // Original rust
//!         a.clone()
//!             .into_iter()
//!             .fold(0f64, |x, y| x + y)
//!             .print(); // 10
//!
//!         // reduce in Peroxide
//!         a.reduce(0f64, |x, y| x + y).print(); // 10
//!     }
//!     ```
//!
//! ### zip_with
//!
//! * `zip_with` is composed of `zip` & `map`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         let b = c!(5,6,7,8);
//!
//!         // Original rust
//!         a.clone()
//!             .into_iter()
//!             .zip(&b)
//!             .map(|(x, y)| x + *y)
//!             .collect::<Vec<f64>>().print();
//!             // [6, 8, 10, 12]
//!
//!         // zip_with in Peroxide
//!         zip_with(|x, y| x + y, &a, &b).print();
//!         // [6, 8, 10, 12]
//!     }
//!     ```
//!
//! ### filter
//!
//! * `filter` is just syntactic sugar for `filter`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         a.filter(|x| x > 2f64).print();
//!         // [3, 4]
//!     }
//!     ```
//!
//! ### take & skip
//!
//! * `take` is syntactic sugar for `take`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         a.take(2).print();
//!         // [1, 2]
//!     }
//!     ```
//!
//! * `skip` is syntactic sugar for `skip`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         a.skip(2).print();
//!         // [3, 4]
//!     }
//!     ```

#[cfg(feature = "openblas")]
extern crate blas;
#[cfg(feature = "openblas")]
use blas::{daxpy, ddot, dnrm2, dscal, idamax};

use operation::extra_ops::Real;
use std::cmp::min;
use std::convert;
use std::f64::MIN;

pub type Vector = Vec<f64>;

/// Functional Programming tools for Vector
pub trait FPVector {
    type Scalar;

    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> Self::Scalar;
    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
        T: convert::Into<Self::Scalar>;
    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar;
    fn filter<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> bool;
    fn take(&self, n: usize) -> Self;
    fn skip(&self, n: usize) -> Self;
}

impl FPVector for Vector {
    type Scalar = f64;

    /// fmap for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// assert_eq!(a.fmap(|x| x*2f64), seq!(2,10,2));
    /// ```
    fn fmap<F>(&self, f: F) -> Vector
    where
        F: Fn(f64) -> f64,
    {
        self.clone().into_iter().map(|x| f(x)).collect::<Vector>()
    }

    /// reduce for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = seq!(1,100,1);
    /// assert_eq!(a.reduce(0, |x,y| x + y), 5050f64);
    /// ```
    fn reduce<F, T>(&self, init: T, f: F) -> f64
    where
        F: Fn(f64, f64) -> f64,
        T: convert::Into<f64>,
    {
        self.clone().into_iter().fold(init.into(), |x, y| f(x, y))
    }

    fn zip_with<F>(&self, f: F, other: &Vector) -> Vector
    where
        F: Fn(f64, f64) -> f64,
    {
        self.into_iter()
            .zip(other)
            .map(|(x, y)| f(*x, *y))
            .collect::<Vector>()
    }

    /// Filter for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// let b = a.filter(|x| x > 3.);
    /// assert_eq!(b, c!(4,5));
    /// ```
    fn filter<F>(&self, f: F) -> Vector
    where
        F: Fn(f64) -> bool,
    {
        self.clone()
            .into_iter()
            .filter(|x| f(*x))
            .collect::<Vector>()
    }

    /// Take for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// let b = a.take(3);
    /// assert_eq!(b, c!(1,2,3));
    /// ```
    fn take(&self, n: usize) -> Vector {
        let mut v = vec![0f64; n];
        for i in 0..n {
            v[i] = self[i];
        }
        return v;
    }

    /// Skip for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// let b = a.skip(3);
    /// assert_eq!(b, c!(4,5));
    /// ```
    fn skip(&self, n: usize) -> Vector {
        let l = self.len();
        let mut v = vec![0f64; l - n];
        for i in n..l {
            v[i - n] = self[i];
        }
        return v;
    }
}

pub fn map<F, T>(f: F, xs: &Vec<T>) -> Vec<T>
where
    F: Fn(T) -> T,
    T: Real + Default,
{
    let l = xs.len();
    let mut result = vec![T::default(); l];
    for i in 0..l {
        result[i] = f(xs[i]);
    }
    result
}

pub fn reduce<F, T>(f: F, init: T, xs: &Vec<T>) -> T
where
    F: Fn(T, T) -> T,
    T: Real,
{
    let mut s = init;
    for i in 0..xs.len() {
        s = f(s, xs[i]);
    }
    s
}

pub fn zip_with<F, T>(f: F, xs: &Vec<T>, ys: &Vec<T>) -> Vec<T>
where
    F: Fn(T, T) -> T,
    T: Real + Default,
{
    let l = min(xs.len(), ys.len());
    let mut result = vec![T::default(); l];
    for i in 0..l {
        result[i] = f(xs[i], ys[i]);
    }
    result
}

/// Some algorithms for Vector
pub trait Algorithm {
    fn rank(&self) -> Vec<usize>;
    fn sign(&self) -> f64;
    fn arg_max(&self) -> usize;
}

impl Algorithm for Vector {
    /// Assign rank
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let v = c!(7, 5, 9, 2, 8);
    /// assert_eq!(v.rank(), vec![2,3,0,4,1]);
    /// ```
    fn rank(&self) -> Vec<usize> {
        let l = self.len();
        let idx = (1..(l + 1)).map(|x| x as usize).collect::<Vec<usize>>();

        let mut vec_tup = self
            .clone()
            .into_iter()
            .zip(idx.clone())
            .collect::<Vec<(f64, usize)>>();
        vec_tup.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap().reverse());
        let indices = vec_tup.into_iter().map(|(_, y)| y).collect::<Vec<usize>>();
        idx.into_iter()
            .map(|x| indices.clone().into_iter().position(|t| t == x).unwrap())
            .collect::<Vec<usize>>()
    }

    /// Sign of Permutation
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,0,2);
    /// let b = c!(1,2,0);
    /// let c = c!(0,1,2);
    ///
    /// assert_eq!((a.sign(), b.sign(), c.sign()), (-1f64, 1f64, 1f64));
    /// ```
    fn sign(&self) -> f64 {
        let l = self.len();
        let mut v = self.clone();
        let mut sgn = 1f64;

        for i in 0..(l - 1) {
            for j in 0..(l - 1 - i) {
                if v[j] > v[j + 1] {
                    sgn *= -1f64;
                    let (a, b) = (v[j], v[j + 1]);
                    v[j] = b;
                    v[j + 1] = a;
                }
            }
        }
        return sgn;
    }

    /// arg max
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let v = c!(1,3,2,4,3,7);
    /// assert_eq!(v.arg_max(),5);
    ///
    /// let v2 = c!(1,3,2,5,6,6);
    /// assert_eq!(v2.arg_max(),4);
    /// ```
    fn arg_max(&self) -> usize {
        match () {
            #[cfg(feature = "native")]
            () => {
                let n_i32 = self.len() as i32;
                let idx: usize;
                unsafe {
                    idx = idamax(n_i32, self, 1) - 1;
                }
                idx
            }
            _ => {
                let v = self.clone();
                let m = self.clone().into_iter().fold(MIN, |x, y| x.max(y));
                v.into_iter().position(|x| x == m).unwrap()
            }
        }
    }
}

/// Convenient Vector Operation trait
pub trait VecOps {
    type Scalar;
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn div(&self, other: &Self) -> Self;
    fn s_add(&self, scala: f64) -> Self;
    fn s_sub(&self, scala: f64) -> Self;
    fn s_mul(&self, scala: f64) -> Self;
    fn s_div(&self, scala: f64) -> Self;
    fn dot(&self, other: &Self) -> Self::Scalar;
    fn norm(&self) -> Self::Scalar;
    fn normalize(&self) -> Self;
}

/// Convenient Vector Operations (No Clone, No Copy)
impl VecOps for Vector {
    type Scalar = f64;

    /// Addition
    fn add(&self, other: &Self) -> Self {
        match () {
            #[cfg(feature = "openblas")]
            () => {
                let n_i32 = self.len() as i32;
                let mut y = other.clone();
                unsafe {
                    daxpy(n_i32, 1f64, self, 1, &mut y, 1);
                }
                y
            }
            _ => self.zip_with(|x, y| x + y, other),
        }
    }

    /// Subtraction
    fn sub(&self, other: &Self) -> Self {
        match () {
            #[cfg(feature = "openblas")]
            () => {
                let n_i32 = self.len() as i32;
                let mut y = self.clone();
                unsafe {
                    daxpy(n_i32, -1f64, other, 1, &mut y, 1);
                }
                y
            }
            _ => self.zip_with(|x, y| x - y, other),
        }
    }

    /// Multiplication
    fn mul(&self, other: &Self) -> Self {
        self.zip_with(|x, y| x * y, other)
    }

    /// Division
    fn div(&self, other: &Self) -> Self {
        self.zip_with(|x, y| x / y, other)
    }

    fn s_add(&self, scala: f64) -> Self {
        match () {
            #[cfg(feature = "openblas")]
            () => {
                let n_i32 = self.len() as i32;
                let mut y = self.clone();
                unsafe {
                    daxpy(n_i32, 1f64, &vec![scala; self.len()], 1, &mut y, 1);
                }
                y
            }
            _ => self.fmap(|x| x + scala),
        }
    }

    fn s_sub(&self, scala: f64) -> Self {
        match () {
            #[cfg(feature = "openblas")]
            () => {
                let n_i32 = self.len() as i32;
                let mut y = self.clone();
                unsafe {
                    daxpy(n_i32, -1f64, &vec![scala; self.len()], 1, &mut y, 1);
                }
                y
            }
            _ => self.fmap(|x| x - scala),
        }
    }

    fn s_mul(&self, scala: f64) -> Self {
        match () {
            #[cfg(feature = "openblas")]
            () => {
                let mut x = self.clone();
                let n_i32 = self.len() as i32;
                unsafe {
                    dscal(n_i32, scala, &mut x, 1);
                }
                x
            }
            _ => self.fmap(|x| scala * x),
        }
    }

    fn s_div(&self, scala: f64) -> Self {
        match () {
            #[cfg(feature = "openblas")]
            () => {
                let mut x = self.clone();
                let n_i32 = self.len() as i32;
                unsafe {
                    dscal(n_i32, 1f64 / scala, &mut x, 1);
                }
                x
            }
            _ => self.fmap(|x| x / scala),
        }
    }

    /// Dot product
    fn dot(&self, other: &Self) -> f64 {
        match () {
            #[cfg(feature = "openblas")]
            () => {
                let n_i32 = self.len() as i32;
                let res: f64;
                unsafe {
                    res = ddot(n_i32, self, 1, other, 1);
                }
                res
            }
            _ => zip_with(|x, y| x * y, &self, other).reduce(0, |x, y| x + y),
        }
    }

    /// Norm
    fn norm(&self) -> f64 {
        match () {
            #[cfg(feature = "openblas")]
            () => {
                let n_i32 = self.len() as i32;
                let res: f64;
                unsafe {
                    res = dnrm2(n_i32, self, 1);
                }
                res
            }
            _ => self.dot(&self).sqrt(),
        }
    }

    /// Normalize
    fn normalize(&self) -> Self {
        let alpha = self.norm();
        self.s_mul(1f64 / alpha)
    }
}

#[cfg(feature = "openblas")]
pub fn blas_daxpy(a: f64, x: &Vec<f64>, y: &mut Vec<f64>) {
    assert_eq!(x.len(), y.len());
    unsafe {
        let n = x.len() as i32;
        daxpy(n, a, x, 1, y, 1)
    }
}

#[cfg(feature = "openblas")]
pub fn blas_daxpy_return(a: f64, x: &Vec<f64>, y: &Vec<f64>) -> Vec<f64> {
    assert_eq!(x.len(), y.len());
    let mut result = y.clone();
    let n = x.len() as i32;
    unsafe {
        daxpy(n, a, x, 1, &mut result, 1);
    }
    result
}
