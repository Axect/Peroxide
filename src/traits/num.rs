//! Missing operations & comprehensive number structures
//!
//! ## `Real` trait
//!
//! * `Real` is a trait for binding `f64`, `AD`
//! * `Real` requires `PowOps, TrigOps, ExpLogOps` & `std::Ops<Self>` & `std::Ops<f64>`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let x_f64 = 2f64;
//!         let x_ad1 = AD1(2f64,1f64);
//!         let x_ad2 = AD2(2f64, 1f64, 0f64);
//!
//!         f(x_f64).print();
//!         f(x_ad1).print();
//!         f(x_ad2).print();
//!     }
//!
//!     fn f<T: Real>(x: T) -> T {
//!         return x.powi(2)
//!     }
//!     ```

use std::ops::{Neg, Add, Sub, Mul, Div};
use crate::structure::ad::AD;

pub trait PowOps: Sized {
    fn powi(&self, n: i32) -> Self;
    fn powf(&self, f: f64) -> Self;
    fn pow(&self, f: Self) -> Self;
    fn sqrt(&self) -> Self {
        self.powf(0.5)
    }
}

pub trait TrigOps: Sized + Div<Output = Self> {
    fn sin_cos(&self) -> (Self, Self);
    fn sinh_cosh(&self) -> (Self, Self);
    fn sin(&self) -> Self {
        let (s, _) = self.sin_cos();
        s
    }
    fn cos(&self) -> Self {
        let (_, c) = self.sin_cos();
        c
    }
    fn tan(&self) -> Self {
        let (s, c) = self.sin_cos();
        s / c
    }
    fn asin(&self) -> Self;
    fn acos(&self) -> Self;
    fn atan(&self) -> Self;
    fn sinh(&self) -> Self {
        let (s, _) = self.sinh_cosh();
        s
    }
    fn cosh(&self) -> Self {
        let (_, c) = self.sinh_cosh();
        c
    }
    fn tanh(&self) -> Self {
        let (s, c) = self.sinh_cosh();
        s / c
    }
    fn asinh(&self) -> Self;
    fn acosh(&self) -> Self;
    fn atanh(&self) -> Self;
}

pub trait ExpLogOps: Sized {
    fn exp(&self) -> Self;
    fn ln(&self) -> Self;
    fn log(&self, base: f64) -> Self;
    fn log2(&self) -> Self {
        self.log(2f64)
    }
    fn log10(&self) -> Self {
        self.log(10f64)
    }
}

pub trait Real:
    PowOps
    + TrigOps
    + ExpLogOps
    + Neg
    + PartialOrd
    + Add<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Sub<Output = Self>
    + Add<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
    + Sub<f64, Output = Self>
    + Clone
    + Copy
{
    fn to_f64(&self) -> f64;
    fn from_f64(f: f64) -> Self;
    fn to_ad(&self) -> AD;
}

// =============================================================================
// Real Traits for f64
// =============================================================================
impl PowOps for f64 {
    fn powi(&self, n: i32) -> Self {
        (*self).powi(n)
    }

    fn powf(&self, f: f64) -> Self {
        (*self).powf(f)
    }

    fn pow(&self, f: f64) -> Self {
        (*self).powf(f)
    }

    fn sqrt(&self) -> Self {
        (*self).sqrt()
    }
}

impl TrigOps for f64 {
    fn sin(&self) -> Self {
        (*self).sin()
    }

    fn cos(&self) -> Self {
        (*self).cos()
    }

    fn tan(&self) -> Self {
        (*self).tan()
    }

    fn asin(&self) -> Self {
        (*self).asin()
    }

    fn acos(&self) -> Self {
        (*self).acos()
    }

    fn atan(&self) -> Self {
        (*self).atan()
    }

    fn sinh(&self) -> Self {
        (*self).sinh()
    }

    fn cosh(&self) -> Self {
        (*self).cosh()
    }

    fn sinh_cosh(&self) -> (Self, Self) {
        ((*self).sinh(), (*self).cosh())
    }

    fn tanh(&self) -> Self {
        (*self).tanh()
    }

    fn asinh(&self) -> Self {
        (*self).asinh()
    }

    fn acosh(&self) -> Self {
        (*self).acosh()
    }

    fn atanh(&self) -> Self {
        (*self).atanh()
    }

    fn sin_cos(&self) -> (Self, Self) {
        (*self).sin_cos()
    }
}

impl ExpLogOps for f64 {
    fn exp(&self) -> Self {
        (*self).exp()
    }

    fn ln(&self) -> Self {
        (*self).ln()
    }

    fn log(&self, base: f64) -> Self {
        (*self).log(base)
    }

    fn log2(&self) -> Self {
        (*self).log2()
    }

    fn log10(&self) -> Self {
        (*self).log10()
    }
}

impl Real for f64 {
    fn to_f64(&self) -> f64 {
        *self
    }

    fn from_f64(f: f64) -> Self {
        f
    }

    fn to_ad(&self) -> AD {
        AD::from(*self)
    }
}

impl Real for AD {
    fn to_f64(&self) -> f64 {
        self.x()
    }

    fn from_f64(f: f64) -> Self {
        AD::from(f)
    }

    fn to_ad(&self) -> AD {
        *self
    }
}
