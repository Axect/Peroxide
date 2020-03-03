//! Missing operations & comprehensive number structures
//!
//! ## `Real` trait
//!
//! * `Real` is a trait for binding `f64`, `Dual`, `HyperDual`
//! * `Real` requires `PowOps, TrigOps, ExpLogOps` & `std::Ops<Self>` & `std::Ops<f64>`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let x_f64 = 2f64;
//!         let x_dual = dual(2,1);
//!         let x_hyper = hyper_dual(2, 1, 0);
//!
//!         f(x_f64).print();
//!         f(x_dual).print();
//!         f(x_hyper).print();
//!     }
//!
//!     fn f<T: Real>(x: T) -> T {
//!         return x.powi(2)
//!     }
//!     ```

use std::ops::{Add, Div, Mul, Sub, Neg};
use structure::dual::Dual;
use structure::hyper_dual::HyperDual;

pub trait PowOps: Sized {
    fn powi(&self, n: i32) -> Self;
    fn powf(&self, f: f64) -> Self;
    fn pow(&self, f: Self) -> Self;
    fn sqrt(&self) -> Self;
}

pub trait TrigOps: Sized {
    fn sin(&self) -> Self;
    fn cos(&self) -> Self;
    fn tan(&self) -> Self;
    fn asin(&self) -> Self;
    fn acos(&self) -> Self;
    fn atan(&self) -> Self;
    fn sinh(&self) -> Self;
    fn cosh(&self) -> Self;
    fn tanh(&self) -> Self;
    fn asinh(&self) -> Self;
    fn acosh(&self) -> Self;
    fn atanh(&self) -> Self;
    fn sin_cos(&self) -> (Self, Self);
}

pub trait ExpLogOps: Sized {
    fn exp(&self) -> Self;
    fn ln(&self) -> Self;
    fn log(&self, base: f64) -> Self;
    fn log2(&self) -> Self;
    fn log10(&self) -> Self;
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
    fn to_dual(&self) -> Dual;
    fn from_dual(d: Dual) -> Self;
    fn to_hyper_dual(&self) -> HyperDual;
    fn from_hyper_dual(h: HyperDual) -> Self;
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

    fn to_dual(&self) -> Dual {
        Dual::new(*self, 0f64)
    }

    fn from_dual(d: Dual) -> Self {
        d.value()
    }

    fn to_hyper_dual(&self) -> HyperDual {
        HyperDual::new(*self, 0f64, 0f64)
    }

    fn from_hyper_dual(h: HyperDual) -> Self {
        h.value()
    }
}

// =============================================================================
// Real trait for Dual
// =============================================================================
impl Real for Dual {
    fn to_f64(&self) -> f64 {
        self.value()
    }

    fn from_f64(f: f64) -> Self {
        Dual::new(f, 0f64)
    }

    fn to_dual(&self) -> Dual {
        *self
    }

    fn from_dual(d: Dual) -> Self {
        d
    }

    fn to_hyper_dual(&self) -> HyperDual {
        HyperDual::new(self.value(), self.slope(), 0f64)
    }

    fn from_hyper_dual(h: HyperDual) -> Self {
        Dual::new(h.value(), h.slope())
    }
}

// =============================================================================
// Real trait for HyperDual
// =============================================================================
impl Real for HyperDual {
    fn to_f64(&self) -> f64 {
        self.value()
    }

    fn from_f64(f: f64) -> Self {
        HyperDual::new(f, 0f64, 0f64)
    }

    fn to_dual(&self) -> Dual {
        Dual::new(self.value(), self.slope())
    }

    fn from_dual(d: Dual) -> Self {
        HyperDual::new(d.value(), d.slope(), 0f64)
    }

    fn to_hyper_dual(&self) -> HyperDual {
        *self
    }

    fn from_hyper_dual(h: HyperDual) -> Self {
        h
    }
}
