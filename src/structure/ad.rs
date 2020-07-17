//! Taylor mode forward automatic differentiation
//!
//! **Caution**: In this documentation, `{i}` means number `i` (ex: AD{2} means AD2)
//!
//! ## Important Features
//!
//! * Can automatic differentiate up to 10th order (`AD1` ~ `AD10`)
//! * All `AD{i}` are in stack (Guarantee high performance)
//! * You can see `AD{i}` via `.print()`
//! * You can use n-th derivative value of `AD{i}` via `.d{i}`
//!
//! ## Implemented Traits for `AD{i}`
//!
//! * `#[derive(Debug, Copy, Clone, PartialEq, Default)]`
//! * `std::fmt::Display`
//! * `From<AD{j}>`, `From<&'a AD{j}>`
//! * `IntoIterator<Item = f64>`
//! * `FromIterator<f64>`
//! * `Index`, `IndexMut`
//! * `std::ops::{Neg, Add, Sub, Mul, Div}`
//! * `peroxide::traits::num::{PowOps, ExpLogOps, TrigOps}`
//!
//! ## Iterator of `AD{i}`
//!
//! There are three iterators.
//!
//! * `ADIntoIter{i}`
//! * `ADIter{i}<'a>`
//! * `ADIterMut{i}<'a>`
//!
//! Each implements `DoubleEndedIterator`, `ExactSizeIterator` also.
//!
//! ## Methods
//!
//! * `new(d0: f64, ... , d{i}: f64) -> AD{i}`
//! * `print(&self)`
//! * `iter(&self) -> ADIter{i}`
//! * `iter_mut(&self) -> ADIterMut{i}`
//! * `len(&self) -> usize`
//!
//! ## Implemented Operations
//!
//! * `Add, Sub, Mul, Div`
//! * `sin, cos, tan`
//! * `sinh, cosh, tanh`
//! * `sin_cos`, `sinh_cosh`
//! * `exp, ln, log, log2, log10`
//! * `powi, powf, sqrt`
//!
//! ## Not yet implemented
//!
//! * `asin`, `acos`, `atan`
//! * `asinh`, `acosh`, `atanh`
//! * `pow`
//!
//! ## Usage
//!
//! ### Construction
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // Declare x where x = 2
//!     let a = AD1::new(2f64, 1f64);
//!     // Declare x^2 where x = 2
//!     let b = AD2::new(4f64, 4f64, 2f64);
//!     // Convert AD1 -> AD2
//!     let c = AD2::from(a);
//!     // Zeros
//!     let d = AD2::default();
//!
//!     assert_eq!(c, AD2::new(2f64, 1f64, 0f64));
//!     assert_eq!(d, AD2::new(0f64, 0f64, 0f64));
//! }
//! ```
//!
//! ### Operation
//!
//! For every binary operation, it returns higher order AD
//! (ex: AD1 + AD2 = AD2)
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = AD1::new(2f64, 1f64);       // x        at x = 2
//!     let b = AD2::new(4f64, 4f64, 2f64); // x^2      at x = 2
//!     let c = a + b;                      // x^2 + x  at x = 2
//!     let d = a * b;                      // x^3      at x = 2
//!     let e = a / b;                      // 1/x      at x = 2
//!     assert_eq!(c, AD2::new(6f64, 5f64, 2f64));
//!     assert_eq!(d, AD2::new(8f64, 12f64, 12f64));
//!     assert_eq!(e, AD2::new(0.5, -0.25, 0.25));
//! }
//! ```
//!

use crate::statistics::ops::C;
use crate::traits::num::{ExpLogOps, PowOps, TrigOps};
use peroxide_ad::{
    ad_display, ad_impl, ad_impl_add, ad_impl_div, ad_impl_double_ended_iter,
    ad_impl_exact_size_iter, ad_impl_explogops, ad_impl_from, ad_impl_from_iter, ad_impl_index,
    ad_impl_into_iter, ad_impl_iter, ad_impl_mul, ad_impl_neg, ad_impl_powops, ad_impl_sub,
    ad_impl_trigops, ad_iter_def, ad_struct_def, ad_impl_ad
};
use std::iter::FromIterator;
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

ad_struct_def!();
ad_display!();
ad_impl!();
ad_impl_from!();
ad_iter_def!();
ad_impl_into_iter!();
ad_impl_iter!();
ad_impl_from_iter!();
ad_impl_double_ended_iter!();
ad_impl_exact_size_iter!();
ad_impl_index!();
ad_impl_neg!();
ad_impl_add!();
ad_impl_sub!();
ad_impl_mul!();
ad_impl_div!();
ad_impl_explogops!();
ad_impl_powops!();
ad_impl_trigops!();
ad_impl_ad!();

pub trait AD:
    std::fmt::Display
    + Clone
    + Copy 
    + PartialEq
    + From<AD1>
    + From<AD2>
    + From<AD3>
    + From<AD4>
    + From<AD5>
    + From<AD6>
    + From<AD7>
    + From<AD8>
    + From<AD9>
    + From<AD10>
    + Into<AD1>
    + Into<AD2>
    + Into<AD3>
    + Into<AD4>
    + Into<AD5>
    + Into<AD6>
    + Into<AD7>
    + Into<AD8>
    + Into<AD9>
    + Into<AD10>
    + IntoIterator<Item=f64>
    + FromIterator<f64>
    + Index<usize>
    + IndexMut<usize>
    + Add<Output=Self>
    + Sub<Output=Self>
    + Mul<Output=Self>
    + Div<Output=Self>
    + PowOps
    + ExpLogOps
    + TrigOps
{
    fn to_ad1(self) -> AD1 {
        self.into()
    }

    fn to_ad2(self) -> AD2 {
        self.into()
    }

    fn to_ad3(self) -> AD3 {
        self.into()
    }

    fn to_ad4(self) -> AD4 {
        self.into()
    }

    fn to_ad5(self) -> AD5 {
        self.into()
    }

    fn to_ad6(self) -> AD6 {
        self.into()
    }

    fn to_ad7(self) -> AD7 {
        self.into()
    }

    fn to_ad8(self) -> AD8 {
        self.into()
    }

    fn to_ad9(self) -> AD9 {
        self.into()
    }

    fn to_ad10(self) -> AD10 {
        self.into()
    }
}
