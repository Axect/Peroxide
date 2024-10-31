//! Pointer toolbox
//!
//! # Redox
//!
//! ## Type
//! ```ignore
//! pub struct Redox<T: Vector> {
//!     data: Box<T>
//! }
//! ```
//!
//! ## Description
//!
//! Operation with `Vec<_>` is too bothered. For example, next code generates error.
//! ```compile_fail
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::prelude::*;
//!
//! fn main() {
//!     let a = c!(1, 2, 3);
//!     assert_eq!(a * 2f64 - 1f64, c!(1, 3, 5));
//! }
//! ```
//!
//! Because we can't implement `Mul<Vec<f64>> for f64` and vice versa.
//! `Redox<T: Vector>` makes the situation easy.
//!
//! ## Usage
//!
//! * `ox()`: `Vector` to `Redox<T: Vector>`
//! * `red()`: `Redox<T: Vector>` to `T` (Ofcourse, `T` should be sized)
//!
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = c!(1, 2, 3);
//!     assert_eq!((a.ox() * 2f64 - 1f64).red(), c!(1, 3, 5));
//! }
//! ```
//!
//! `ox()` and `red()` come from oxidation and reduction.
use crate::structure::ad::AD;
use crate::structure::matrix::{Matrix, Shape};
use crate::structure::sparse::SPMatrix;
use crate::traits::{
    fp::FPVector,
    math::{LinearOp, Vector},
    matrix::MatrixTrait,
};
use std::ops::{Add, Deref, Div, Mul, Sub};

// =============================================================================
// Redox Structure
// =============================================================================
#[derive(Debug)]
pub struct Redox<T: Vector> {
    data: Box<T>,
}

impl<T: Vector> Deref for Redox<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

pub trait RedoxCommon {
    type ToRedox;
    fn from_vec(vec: Self::ToRedox) -> Self;
    fn red(self) -> Self::ToRedox;
}

impl RedoxCommon for Redox<Vec<f64>> {
    type ToRedox = Vec<f64>;
    fn from_vec(vec: Self::ToRedox) -> Self {
        Self {
            data: Box::new(vec),
        }
    }

    fn red(self) -> Self::ToRedox {
        (*self).to_vec()
    }
}

impl RedoxCommon for Redox<Vec<AD>> {
    type ToRedox = Vec<AD>;
    fn from_vec(vec: Self::ToRedox) -> Self {
        Self {
            data: Box::new(vec),
        }
    }

    fn red(self) -> Self::ToRedox {
        (*self).to_vec()
    }
}

// =============================================================================
// Oxide trait
// =============================================================================
pub trait Oxide: Vector {
    fn ox(self) -> Redox<Self>
    where
        Self: Sized;
}

// =============================================================================
// Reference Arithmetic
// =============================================================================
impl<T: Vector> Add<Redox<T>> for Redox<T> {
    type Output = Self;

    fn add(self, rhs: Redox<T>) -> Self::Output {
        Redox {
            data: Box::new(self.add_vec(&rhs.data)),
        }
    }
}

impl<T: Vector + FPVector> Sub<Redox<T>> for Redox<T>
where
    <T as FPVector>::Scalar: Sub<Output = <T as FPVector>::Scalar>,
{
    type Output = Self;

    fn sub(self, rhs: Redox<T>) -> Self::Output {
        Redox {
            data: Box::new(self.zip_with(|x, y| x - y, &rhs.data)),
        }
    }
}

impl<T: Vector + FPVector> Mul<Redox<T>> for Redox<T>
where
    <T as FPVector>::Scalar: Mul<Output = <T as FPVector>::Scalar>,
{
    type Output = Self;

    fn mul(self, rhs: Redox<T>) -> Self::Output {
        Redox {
            data: Box::new(self.zip_with(|x, y| x * y, &rhs.data)),
        }
    }
}

impl<T: Vector + FPVector> Div<Redox<T>> for Redox<T>
where
    <T as FPVector>::Scalar: Div<Output = <T as FPVector>::Scalar>,
{
    type Output = Self;

    fn div(self, rhs: Redox<T>) -> Self::Output {
        Redox {
            data: Box::new(self.zip_with(|x, y| x / y, &rhs.data)),
        }
    }
}

impl<T: Vector + FPVector> Add<f64> for Redox<T>
where
    <T as FPVector>::Scalar: Add<f64, Output = <T as FPVector>::Scalar>,
{
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        Redox {
            data: Box::new(self.fmap(|x| x + rhs)),
        }
    }
}

impl<T: Vector + FPVector> Sub<f64> for Redox<T>
where
    <T as FPVector>::Scalar: Sub<f64, Output = <T as FPVector>::Scalar>,
{
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        Redox {
            data: Box::new(self.fmap(|x| x - rhs)),
        }
    }
}

impl<T: Vector + FPVector> Mul<f64> for Redox<T>
where
    <T as FPVector>::Scalar: Mul<f64, Output = <T as FPVector>::Scalar>,
{
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Redox {
            data: Box::new(self.fmap(|x| x * rhs)),
        }
    }
}

impl<T: Vector + FPVector> Div<f64> for Redox<T>
where
    <T as FPVector>::Scalar: Div<f64, Output = <T as FPVector>::Scalar>,
{
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Redox {
            data: Box::new(self.fmap(|x| x / rhs)),
        }
    }
}

impl Mul<Redox<Vec<f64>>> for Matrix {
    type Output = Redox<Vec<f64>>;

    fn mul(self, rhs: Redox<Vec<f64>>) -> Self::Output {
        Redox {
            data: Box::new(self.apply(&*rhs)),
        }
    }
}

impl Mul<Redox<Vec<f64>>> for &Matrix {
    type Output = Redox<Vec<f64>>;

    fn mul(self, rhs: Redox<Vec<f64>>) -> Self::Output {
        Redox {
            data: Box::new(self.apply(&*rhs)),
        }
    }
}

/// Matrix multiplication with Redox
impl Mul<Redox<Vec<f64>>> for SPMatrix {
    type Output = Redox<Vec<f64>>;
    fn mul(self, rhs: Redox<Vec<f64>>) -> Self::Output {
        Redox {
            data: Box::new(self.apply(&rhs.data)),
        }
    }
}

impl Mul<Redox<Vec<f64>>> for &SPMatrix {
    type Output = Redox<Vec<f64>>;

    fn mul(self, rhs: Redox<Vec<f64>>) -> Self::Output {
        Redox {
            data: Box::new(self.apply(&rhs.data)),
        }
    }
}

// =============================================================================
// Pointer for Matrix
// =============================================================================
/// Pointer for col or row
pub trait MatrixPtr {
    unsafe fn row_ptr(&self, idx: usize) -> Vec<*const f64>;
    unsafe fn col_ptr(&self, idx: usize) -> Vec<*const f64>;
}

impl MatrixPtr for Matrix {
    /// Row pointer
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = ml_matrix("1 2;3 4");
    ///     unsafe {
    ///         let b = a.row_ptr(1);
    ///         let b_vec = ptr_to_vec(&b);
    ///         assert_eq!(b_vec, vec![3f64, 4f64]);
    ///     }
    /// }
    /// ```
    unsafe fn row_ptr(&self, idx: usize) -> Vec<*const f64> {
        assert!(idx < self.col, "Index out of range");
        match self.shape {
            Shape::Row => {
                let mut v: Vec<*const f64> = vec![&0f64; self.col];
                let start_idx = idx * self.col;
                let p = self.ptr();
                for (i, j) in (start_idx..start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Col => {
                let mut v: Vec<*const f64> = vec![&0f64; self.col];
                let p = self.ptr();
                for (i, elem) in v.iter_mut().enumerate() {
                    *elem = p.add(idx + i * self.row);
                }
                v
            }
        }
    }

    unsafe fn col_ptr(&self, idx: usize) -> Vec<*const f64> {
        assert!(idx < self.col, "Index out of range");
        match self.shape {
            Shape::Col => {
                let mut v: Vec<*const f64> = vec![&0f64; self.row];
                let start_idx = idx * self.row;
                let p = self.ptr();
                for (i, j) in (start_idx..start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Row => {
                let mut v: Vec<*const f64> = vec![&0f64; self.row];
                let p = self.ptr();
                for (i, elem) in v.iter_mut().enumerate() {
                    *elem = p.add(idx + i * self.col);
                }
                v
            }
        }
    }
}
