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

use crate::structure::ad::AD;
use peroxide_num::{ExpLogOps, PowOps, TrigOps};
use std::ops::{Add, Div, Mul, Neg, Sub};

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
