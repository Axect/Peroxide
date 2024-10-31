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
//!     use peroxide::fuga::*;
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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!     }
//!     ```
//!
//! ## From ranges to `Vec<f64>`
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
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = seq(1, 4, 1);
//!         a.print();
//!         // [1, 2, 3, 4]
//!     }
//!     ```
//!
//! ## `Vec<f64>` Operation
//!
//! There are two ways to do vector operations.
//!
//! * Use functional programming tools
//! * Use redox
//!
//! Here, I explain second method - for functional programming, see below.
//!
//! To use redox, you only need to understand two things - `ox()`, `red()`.
//!
//! * `ox()` : Makes vector to `Redox<T: Vector>`
//! * `red()` : Makes `Redox<T: Vector>` to vector.
//!
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = c!(1, 2, 3);
//!     assert_eq!((a.ox() * 2f64 - 1f64).red(), c!(1f64, 3f64, 5f64));
//! }
//! ```
//!
//! ## Concatenation
//!
//! There are two concatenation operations.
//!
//! * `cat(T, Vec<T>) -> Vec<f64>`
//! * `concat(Vec<T>, Vec<T>) -> Vec<T>`
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         cat(0f64, &a).print();
//!         // [0, 1, 2, 3, 4]
//!
//!         let b = c!(5,6,7,8);
//!         concat(&a, &b).print();
//!         // [1, 2, 3, 4, 5, 6, 7, 8]
//!     }
//!     ```
//!
//! ## Conversion to Matrix
//!
//! There are two ways to convert `Vec<f64>` to `Matrix`.
//!
//! * `into(self) -> Matrix` : `Vec<f64>` to column matrix
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         let a_col: Matrix = a.into();
//!         let m_col = matrix(c!(1,2,3,4), 4, 1, Col); // (4,1) Matrix
//!         assert_eq!(a_col, m_col);
//!
//!         let m_row = matrix(c!(1,2,3,4), 1, 4, Row); // (1,4) Matrix
//!         assert_eq!(a_col.t(), m_row);
//!     }
//!     ```
//!
//! # Functional Programming {#functional}
//!
//! ## FP for `Vec<f64>`
//!
//! * There are some functional programming tools for `Vec<f64>`
//!
//! ### fmap
//!
//! * `fmap` is syntactic sugar for `map`
//! * But different to original `map` - Only `f64 -> f64` allowed.
//!
//!     ```rust
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let a = c!(1,2,3,4);
//!         a.skip(2).print();
//!         // [3, 4]
//!     }
//!     ```

#[cfg(feature = "O3")]
extern crate blas;
#[cfg(feature = "parallel")]
use crate::rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    IntoParallelRefMutIterator, ParallelIterator,
};
#[cfg(feature = "O3")]
use blas::{daxpy, ddot, dnrm2, idamax};

use crate::structure::matrix::{matrix, Matrix, Row};
use crate::traits::{
    fp::FPVector,
    general::Algorithm,
    math::{InnerProduct, LinearOp, Norm, Normed, Vector, VectorProduct},
    mutable::MutFP,
    pointer::{Oxide, Redox, RedoxCommon},
};
#[cfg(feature = "parallel")]
use crate::traits::{
    fp::ParallelFPVector,
    math::{ParallelInnerProduct, ParallelNormed, ParallelVectorProduct},
    mutable::ParallelMutFP,
};
use std::cmp::min;

impl FPVector for Vec<f64> {
    type Scalar = f64;

    /// fmap for `Vec<f64>`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     assert_eq!(a.fmap(|x| x*2f64), seq!(2,10,2));
    /// }
    /// ```
    fn fmap<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
    {
        let mut v = self.clone();
        v.iter_mut().for_each(|x| *x = f(*x));
        v
    }

    /// reduce for `Vec<f64>`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = seq!(1,100,1);
    ///     assert_eq!(a.reduce(0, |x,y| x + y), 5050f64);
    /// }
    /// ```
    fn reduce<F, T>(&self, init: T, f: F) -> f64
    where
        F: Fn(f64, f64) -> f64,
        T: Into<f64> + Copy,
    {
        self.iter().fold(init.into(), |x, &y| f(x, y))
    }

    fn zip_with<F>(&self, f: F, other: &Vec<f64>) -> Vec<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        self.iter()
            .zip(other)
            .map(|(x, y)| f(*x, *y))
            .collect::<Vec<f64>>()
    }

    /// Filter for `Vec<f64>`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     let b = a.filter(|x| x > 3.);
    ///     assert_eq!(b, c!(4,5));
    /// }
    /// ```
    fn filter<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(f64) -> bool,
    {
        self.clone()
            .into_iter()
            .filter(|x| f(*x))
            .collect::<Vec<f64>>()
    }

    /// Take for `Vec<f64>`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     let b = a.take(3);
    ///     assert_eq!(b, c!(1,2,3));
    /// }
    /// ```
    fn take(&self, n: usize) -> Vec<f64> {
        let mut v = vec![0f64; n];
        v[..n].copy_from_slice(&self[..n]);
        v
    }

    /// Skip for `Vec<f64>`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     let b = a.skip(3);
    ///     assert_eq!(b, c!(4,5));
    /// }
    /// ```
    fn skip(&self, n: usize) -> Vec<f64> {
        let l = self.len();
        let mut v = vec![0f64; l - n];
        for (i, j) in (n..l).enumerate() {
            v[i] = self[j];
        }
        v
    }

    fn sum(&self) -> f64 {
        self.iter().sum()
    }

    fn prod(&self) -> f64 {
        self.iter().product()
    }
}

#[cfg(feature = "parallel")]
impl ParallelFPVector for Vec<f64> {
    type Scalar = f64;

    /// par_fmap for `Vec<f64>`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    /// use peroxide::traits::fp::ParallelFPVector;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     assert_eq!(a.par_fmap(|x| x*2f64), seq!(2,10,2));
    /// }
    /// ```
    fn par_fmap<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(f64) -> f64 + Sync + Send,
    {
        let mut v = self.clone();
        v.par_iter_mut().for_each(|x| *x = f(*x));
        v
    }

    /// reduce for `Vec<f64>`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    /// use peroxide::traits::fp::ParallelFPVector;
    ///
    /// fn main() {
    ///     let a = seq!(1,100,1);
    ///     assert_eq!(a.par_reduce(0, |x,y| x + y), 5050f64);
    /// }
    /// ```
    fn par_reduce<F, T>(&self, _init: T, f: F) -> f64
    where
        F: Fn(f64, f64) -> f64 + Send + Sync,
        T: Into<f64> + Send + Sync + Copy,
    {
        self.par_iter()
            .cloned()
            .reduce_with(|x, y| f(x, y))
            .expect("Unable to perform parallel reduce operation")
    }

    fn par_zip_with<F>(&self, f: F, other: &Vec<f64>) -> Vec<f64>
    where
        F: Fn(f64, f64) -> f64 + Send + Sync,
    {
        self.par_iter()
            .zip(other)
            .map(|(x, y)| f(*x, *y))
            .collect::<Vec<f64>>()
    }

    /// Filter for `Vec<f64>`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    /// use peroxide::traits::fp::ParallelFPVector;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5);
    ///     let b = a.par_filter(|x| x > 3.);
    ///     assert_eq!(b, c!(4,5));
    /// }
    /// ```
    fn par_filter<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(f64) -> bool + Send + Sync,
    {
        self.clone()
            .into_par_iter()
            .filter(|x| f(*x))
            .collect::<Vec<f64>>()
    }
}

/// Explicit version of `map`
pub fn map<F, T>(f: F, xs: &[T]) -> Vec<T>
where
    F: Fn(T) -> T,
    T: Default + Copy,
{
    let l = xs.len();
    let mut result = vec![T::default(); l];
    for i in 0..l {
        result[i] = f(xs[i]);
    }
    result
}

/// Explicit version of `map` in parallel
#[cfg(feature = "parallel")]
pub fn par_map<F, T>(f: F, xs: &[T]) -> Vec<T>
where
    F: Fn(T) -> T + Send + Sync,
    T: Default + Copy + Send + Sync,
{
    let mut result = vec![T::default(); xs.len()];
    result.par_iter_mut().for_each(|x| *x = f(*x));
    result
}

/// Explicit version of `reduce`
pub fn reduce<F, T>(f: F, init: T, xs: &[T]) -> T
where
    F: Fn(T, T) -> T,
    T: Copy,
{
    let mut s = init;
    for &x in xs {
        s = f(s, x);
    }
    s
}

/// Explicit version of `zip_with`
pub fn zip_with<F, T>(f: F, xs: &[T], ys: &[T]) -> Vec<T>
where
    F: Fn(T, T) -> T,
    T: Default + Copy,
{
    let l = min(xs.len(), ys.len());
    let mut result = vec![T::default(); l];
    for i in 0..l {
        result[i] = f(xs[i], ys[i]);
    }
    result
}

/// Explicit version of `zip_with` in parallel
#[cfg(feature = "parallel")]
pub fn par_zip_with<F, T>(f: F, xs: &[T], ys: &[T]) -> Vec<T>
where
    F: Fn(T, T) -> T + Send + Sync,
    T: Default + Copy + Send + Sync,
{
    xs.par_iter()
        .zip(ys.into_par_iter())
        .map(|(x, y)| f(*x, *y))
        .collect::<Vec<T>>()
}

impl MutFP for Vec<f64> {
    type Scalar = f64;

    fn mut_map<F>(&mut self, f: F)
    where
        F: Fn(Self::Scalar) -> Self::Scalar,
    {
        for x in self.iter_mut() {
            *x = f(*x);
        }
    }

    fn mut_zip_with<F>(&mut self, f: F, other: &Self)
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
    {
        for i in 0..self.len() {
            self[i] = f(self[i], other[i]);
        }
    }
}

#[cfg(feature = "parallel")]
impl ParallelMutFP for Vec<f64> {
    type Scalar = f64;

    fn par_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Self::Scalar) -> Self::Scalar + Send + Sync,
    {
        self.par_iter_mut().for_each(|x| *x = f(*x));
    }

    fn par_mut_zip_with<F>(&mut self, f: F, other: &Self)
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar + Send + Sync,
    {
        self.par_iter_mut()
            .zip(other.into_par_iter())
            .for_each(|(x, y)| *x = f(*x, *y));
    }
}

impl Algorithm for Vec<f64> {
    type Scalar = f64;

    /// Assign rank
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let v = c!(7, 5, 9, 2, 8);
    ///     assert_eq!(v.rank(), vec![2,3,0,4,1]);
    /// }
    /// ```
    fn rank(&self) -> Vec<usize> {
        let l = self.len();
        let idx = (1..(l + 1)).collect::<Vec<usize>>();

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
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,0,2);
    ///     let b = c!(1,2,0);
    ///     let c = c!(0,1,2);
    ///
    ///     assert_eq!((a.sign(), b.sign(), c.sign()), (-1f64, 1f64, 1f64));
    /// }
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
        sgn
    }

    /// arg max
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let v = c!(1,3,2,4,3,7);
    ///     assert_eq!(v.arg_max(),5);
    ///
    ///     let v2 = c!(1,3,2,5,6,6);
    ///     assert_eq!(v2.arg_max(),4);
    /// }
    /// ```
    fn arg_max(&self) -> usize {
        #[cfg(feature = "O3")]
        {
            let n_i32 = self.len() as i32;
            let idx: usize;
            unsafe {
                idx = idamax(n_i32, self, 1) - 1;
            }
            idx
        }
        #[cfg(not(feature = "O3"))]
        {
            self.into_iter()
                .enumerate()
                .fold((0usize, f64::MIN), |acc, (ics, &val)| {
                    if acc.1 < val {
                        (ics, val)
                    } else {
                        acc
                    }
                })
                .0
        }
    }

    /// arg min
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let v = c!(1,3,2,4,3,7);
    ///     assert_eq!(v.arg_min(),0);
    ///
    ///     let v2 = c!(4,3,2,5,1,6);
    ///     assert_eq!(v2.arg_min(),4);
    /// }
    fn arg_min(&self) -> usize {
        self.iter()
            .enumerate()
            .fold((0usize, f64::MAX), |acc, (ics, &val)| {
                if acc.1 > val {
                    (ics, val)
                } else {
                    acc
                }
            })
            .0
    }

    fn max(&self) -> f64 {
        #[cfg(feature = "O3")]
        {
            let n_i32 = self.len() as i32;
            let idx: usize;
            unsafe {
                idx = idamax(n_i32, self, 1) - 1;
            }
            self[idx]
        }
        #[cfg(not(feature = "O3"))]
        {
            self.into_iter()
                .fold(f64::MIN, |acc, &val| if acc < val { val } else { acc })
        }
    }

    fn min(&self) -> f64 {
        self.iter()
            .fold(f64::MAX, |acc, &val| if acc > val { val } else { acc })
    }

    fn swap_with_perm(&mut self, p: &Vec<(usize, usize)>) {
        for (i, j) in p.iter() {
            self.swap(*i, *j);
        }
    }
}

impl Vector for Vec<f64> {
    type Scalar = f64;

    fn add_vec(&self, rhs: &Self) -> Self {
        self.zip_with(|x, y| x + y, rhs)
    }

    fn sub_vec(&self, rhs: &Self) -> Self {
        self.zip_with(|x, y| x - y, rhs)
    }

    fn mul_scalar(&self, rhs: Self::Scalar) -> Self {
        let alpha: f64 = rhs;
        self.fmap(|x| x * alpha)
    }
}

impl Normed for Vec<f64> {
    type UnsignedScalar = f64;
    fn norm(&self, kind: Norm) -> f64 {
        match kind {
            Norm::L1 => self.iter().map(|x| x.abs()).sum(),
            Norm::L2 => {
                #[cfg(feature = "O3")]
                {
                    let n_i32 = self.len() as i32;
                    let res: f64;
                    unsafe {
                        res = dnrm2(n_i32, self, 1);
                    }
                    res
                }
                #[cfg(not(feature = "O3"))]
                {
                    self.iter().map(|x| x.powi(2)).sum::<f64>().sqrt()
                }
            }
            Norm::Lp(p) => {
                assert!(
                    p >= 1f64,
                    "lp norm is only defined for p>=1, the given value was p={}",
                    p
                );
                self.iter().map(|x| x.powf(p)).sum::<f64>().powf(1f64 / p)
            }
            Norm::LInf => self.iter().fold(0f64, |x, y| x.max(y.abs())),
            Norm::F => unimplemented!(),
            Norm::Lpq(_, _) => unimplemented!(),
        }
    }
    fn normalize(&self, kind: Norm) -> Self
    where
        Self: Sized,
    {
        let denom = self.norm(kind);
        self.fmap(|x| x / denom)
    }
}

#[cfg(feature = "parallel")]
impl ParallelNormed for Vec<f64> {
    type UnsignedScalar = f64;

    fn par_norm(&self, kind: Norm) -> f64 {
        match kind {
            Norm::L1 => self.iter().map(|x| x.abs()).sum(),
            Norm::L2 => {
                #[cfg(feature = "O3")]
                {
                    let n_i32 = self.len() as i32;
                    let res: f64;
                    unsafe {
                        res = dnrm2(n_i32, self, 1);
                    }
                    res
                }
                #[cfg(not(feature = "O3"))]
                {
                    self.par_iter()
                        .map(|x| x.powi(2))
                        .sum::<Self::UnsignedScalar>()
                        .sqrt()
                }
            }
            Norm::Lp(p) => {
                assert!(
                    p >= 1f64,
                    "lp norm is only defined for p>=1, the given value was p={}",
                    p
                );
                self.par_iter()
                    .map(|x| x.powf(p))
                    .sum::<Self::UnsignedScalar>()
                    .powf(1_f64 / p)
            }
            Norm::LInf => self
                .par_iter()
                .fold(|| 0f64, |x, y| x.max(y.abs()))
                .reduce(|| 0f64, |acc1, acc2| acc1.max(acc2.abs())),
            Norm::F => unimplemented!(),
            Norm::Lpq(_, _) => unimplemented!(),
        }
    }
}

impl InnerProduct for Vec<f64> {
    fn dot(&self, rhs: &Self) -> f64 {
        #[cfg(feature = "O3")]
        {
            let n_i32 = self.len() as i32;
            let res: f64;
            unsafe {
                res = ddot(n_i32, self, 1, rhs, 1);
            }
            res
        }
        #[cfg(not(feature = "O3"))]
        {
            self.iter()
                .zip(rhs.iter())
                .fold(0f64, |x, (y1, y2)| x + y1 * y2)
        }
    }
}

#[cfg(feature = "parallel")]
impl ParallelInnerProduct for Vec<f64> {
    fn par_dot(&self, rhs: &Self) -> f64 {
        #[cfg(feature = "O3")]
        {
            let n_i32 = self.len() as i32;
            let res: f64;
            unsafe {
                res = ddot(n_i32, self, 1, rhs, 1);
            }
            res
        }
        #[cfg(not(feature = "O3"))]
        {
            self.par_iter()
                .zip(rhs.into_par_iter())
                .fold(|| 0_f64, |x, (y1, y2)| x + y1 * y2)
                .sum::<f64>()
        }
    }
}

impl LinearOp<Vec<f64>, f64> for Vec<f64> {
    fn apply(&self, rhs: &Vec<f64>) -> f64 {
        self.dot(rhs)
    }
}

impl Oxide for Vec<f64> {
    fn ox(self) -> Redox<Vec<f64>> {
        Redox::<Vec<f64>>::from_vec(self)
    }
}

impl VectorProduct for Vec<f64> {
    fn cross(&self, other: &Self) -> Self {
        assert_eq!(self.len(), other.len());
        // 2D cross product is ill-defined
        match self.len() {
            2 => {
                let mut v = vec![0f64; 1];
                v[0] = self[0] * other[1] - self[1] * other[0];
                v
            }
            3 => {
                let mut v = vec![0f64; 3];
                v[0] = self[1] * other[2] - self[2] * other[1];
                v[1] = self[2] * other[0] - self[0] * other[2];
                v[2] = self[0] * other[1] - self[1] * other[0];
                v
            }
            _ => {
                panic!("Cross product can be defined only in 2 or 3 dimension")
            }
        }
    }

    fn outer(&self, other: &Self) -> Matrix {
        let m: Matrix = self.into();
        let n = matrix(other.to_owned(), 1, other.len(), Row);
        m * n
    }
}

#[cfg(feature = "parallel")]
impl ParallelVectorProduct for Vec<f64> {
    fn par_cross(&self, other: &Self) -> Self {
        assert_eq!(self.len(), other.len());
        // 2D cross product is ill-defined
        match self.len() {
            2 => {
                let mut v = vec![0f64; 1];
                v[0] = self[0] * other[1] - self[1] * other[0];
                v
            }
            3 => {
                let v = (0..3)
                    .into_par_iter()
                    .map(|index| {
                        self[(index + 1) % 3] * other[(index + 2) % 3]
                            - self[(index + 2) % 3] * other[(index + 1) % 3]
                    })
                    .collect::<Vec<f64>>();
                v
            }
            _ => {
                panic!("Cross product can be defined only in 2 or 3 dimension")
            }
        }
    }
}

// /// Convenient Vec<f64> Operations (No Clone, No Copy)
// impl VecOps for Vec<f64> {
//     fn s_add(&self, scala: f64) -> Self {
//         match () {
//             #[cfg(feature = "O3")]
//             () => {
//                 match self.len() {
//                     n if n % 8 == 0 => {
//                         let mut z = vec![0f64; n];
//                         self.chunks_exact(8)
//                             .map(f64x8::from_slice_unaligned)
//                             .map(|x| x + scala)
//                             .for_each(|x| x.write_to_slice_unaligned(&mut z));
//                         z
//                     }
//                     n if n % 4 == 0 => {
//                         let mut z = vec![0f64; n];
//                         self.chunks_exact(4)
//                             .map(f64x4::from_slice_unaligned)
//                             .map(|x| x + scala)
//                             .for_each(|x| x.write_to_slice_unaligned(&mut z));
//                         z
//                     }
//                     _ => self.fmap(|x| x + scala)
//                 }
//             }
//             _ => self.fmap(|x| x + scala),
//         }
//     }
//
//     fn s_sub(&self, scala: f64) -> Self {
//         match () {
//             #[cfg(feature = "O3")]
//             () => {
//                 match self.len() {
//                     n if n % 8 == 0 => {
//                         let mut z = vec![0f64; n];
//                         self.chunks_exact(8)
//                             .map(f64x8::from_slice_unaligned)
//                             .map(|x| x - scala)
//                             .for_each(|x| x.write_to_slice_unaligned(&mut z));
//                         z
//                     }
//                     n if n % 4 == 0 => {
//                         let mut z = vec![0f64; n];
//                         self.chunks_exact(4)
//                             .map(f64x4::from_slice_unaligned)
//                             .map(|x| x - scala)
//                             .for_each(|x| x.write_to_slice_unaligned(&mut z));
//                         z
//                     }
//                     _ => self.fmap(|x| x - scala)
//                 }
//             }
//             _ => self.fmap(|x| x - scala),
//         }
//     }
//
//     fn s_mul(&self, scala: f64) -> Self {
//         match () {
//             #[cfg(feature = "O3")]
//             () => {
//                 match self.len() {
//                     n if n % 8 == 0 => {
//                         let mut z = vec![0f64; n];
//                         self.chunks_exact(8)
//                             .map(f64x8::from_slice_unaligned)
//                             .map(|x| x * scala)
//                             .for_each(|x| x.write_to_slice_unaligned(&mut z));
//                         z
//                     }
//                     n if n % 4 == 0 => {
//                         let mut z = vec![0f64; n];
//                         self.chunks_exact(4)
//                             .map(f64x4::from_slice_unaligned)
//                             .map(|x| x * scala)
//                             .for_each(|x| x.write_to_slice_unaligned(&mut z));
//                         z
//                     }
//                     _ => self.fmap(|x| x * scala)
//                 }
//             }
//             _ => self.fmap(|x| scala * x),
//         }
//     }
//
//     fn s_div(&self, scala: f64) -> Self {
//         match () {
//             #[cfg(feature = "O3")]
//             () => {
//                 match self.len() {
//                     n if n % 8 == 0 => {
//                         let mut z = vec![0f64; n];
//                         self.chunks_exact(8)
//                             .map(f64x8::from_slice_unaligned)
//                             .map(|x| x / scala)
//                             .for_each(|x| x.write_to_slice_unaligned(&mut z));
//                         z
//                     }
//                     n if n % 4 == 0 => {
//                         let mut z = vec![0f64; n];
//                         self.chunks_exact(4)
//                             .map(f64x4::from_slice_unaligned)
//                             .map(|x| x / scala)
//                             .for_each(|x| x.write_to_slice_unaligned(&mut z));
//                         z
//                     }
//                     _ => self.fmap(|x| x / scala)
//                 }
//             }
//             _ => self.fmap(|x| x / scala),
//         }
//     }
//
//     /// Dot product
//     fn dot(&self, other: &Self) -> f64 {
//         match () {
//             #[cfg(feature = "O3")]
//             () => {
//                 let n_i32 = self.len() as i32;
//                 let res: f64;
//                 unsafe {
//                     res = ddot(n_i32, self, 1, other, 1);
//                 }
//                 res
//             }
//             _ => zip_with(|x, y| x * y, &self, other).reduce(0, |x, y| x + y),
//         }
//     }
//
//     fn sum(&self) -> Self::Scalar {
//         match () {
//             #[cfg(feature = "O3")]
//             () => {
//                 let n_i32 = self.len() as i32;
//                 let res: f64;
//                 unsafe {
//                     res = dasum(n_i32, self, 1);
//                 }
//                 res
//             }
//             _ => reduce(|x, y| x + y, 0f64, self)
//         }
//     }

#[cfg(feature = "O3")]
pub fn blas_daxpy(a: f64, x: &[f64], y: &mut [f64]) {
    assert_eq!(x.len(), y.len());
    unsafe {
        let n = x.len() as i32;
        daxpy(n, a, x, 1, y, 1)
    }
}

#[cfg(feature = "O3")]
pub fn blas_daxpy_return(a: f64, x: &[f64], y: &[f64]) -> Vec<f64> {
    assert_eq!(x.len(), y.len());
    let mut result = y.to_vec();
    let n = x.len() as i32;
    unsafe {
        daxpy(n, a, x, 1, &mut result, 1);
    }
    result
}
