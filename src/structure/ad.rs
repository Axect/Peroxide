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
//! ### Generic
//!
//! * All of `AD{i}` implements `AD` trait
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = AD1::new(2f64, 1f64);
//!     let b = AD2::new(4f64, 4f64, 2f64);
//!     assert_eq!(f(a, b), AD1::new(6f64, 5f64));
//! }
//!
//! fn f<T: AD, S: AD>(x: T, y: S) -> T {
//!     T::from(x.to_ad2() + y.to_ad2())
//! }
//! ```

use crate::statistics::ops::C;
use crate::traits::{
    num::{ExpLogOps, PowOps, TrigOps},
    stable::StableFn,
};
use peroxide_ad::{
    ad_display, ad_impl, ad_impl_ad, ad_impl_add, ad_impl_div, ad_impl_double_ended_iter,
    ad_impl_exact_size_iter, ad_impl_explogops, ad_impl_from, ad_impl_from_iter, ad_impl_index,
    ad_impl_into_iter, ad_impl_iter, ad_impl_mul, ad_impl_neg, ad_impl_powops, ad_impl_sub,
    ad_impl_trigops, ad_iter_def, ad_struct_def, ad_impl_from_type, ad_impl_add_f64, ad_impl_sub_f64,
    ad_impl_mul_f64, ad_impl_div_f64, f64_impl_add_ad, f64_impl_sub_ad, f64_impl_mul_ad, f64_impl_div_ad,
    f64_impl_from_ad, ad_impl_stable_fn, def_ad
};
use std::iter::FromIterator;
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};
use std::marker::PhantomData;

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
def_ad!();
ad_impl_from_type!(f64);
// ad_impl_from_type!(f32);
// ad_impl_from_type!(i64);
// ad_impl_from_type!(i32);
ad_impl_add_f64!();
ad_impl_sub_f64!();
ad_impl_mul_f64!();
ad_impl_div_f64!();
f64_impl_add_ad!();
f64_impl_sub_ad!();
f64_impl_mul_ad!();
f64_impl_div_ad!();
f64_impl_from_ad!();
ad_impl_stable_fn!();

// pub trait AD:
//     std::fmt::Display
//     + Clone
//     + Copy
//     + PartialEq
//     + From<AD1>
//     + From<AD2>
//     + From<AD3>
//     + From<AD4>
//     + From<AD5>
//     + From<AD6>
//     + From<AD7>
//     + From<AD8>
//     + From<AD9>
//     + From<AD10>
//     + Into<AD1>
//     + Into<AD2>
//     + Into<AD3>
//     + Into<AD4>
//     + Into<AD5>
//     + Into<AD6>
//     + Into<AD7>
//     + Into<AD8>
//     + Into<AD9>
//     + Into<AD10>
//     + From<f64>
//     + Into<f64>
//     + Add<Output = Self>
//     + Sub<Output = Self>
//     + Mul<Output = Self>
//     + Div<Output = Self>
//     + Add<f64, Output = Self>
//     + Sub<f64, Output = Self>
//     + Mul<f64, Output = Self>
//     + Div<f64, Output = Self>
//     + PowOps
//     + ExpLogOps
//     + TrigOps
// {
//     fn to_ad1(self) -> AD1 {
//         self.into()
//     }
//
//     fn to_ad2(self) -> AD2 {
//         self.into()
//     }
//
//     fn to_ad3(self) -> AD3 {
//         self.into()
//     }
//
//     fn to_ad4(self) -> AD4 {
//         self.into()
//     }
//
//     fn to_ad5(self) -> AD5 {
//         self.into()
//     }
//
//     fn to_ad6(self) -> AD6 {
//         self.into()
//     }
//
//     fn to_ad7(self) -> AD7 {
//         self.into()
//     }
//
//     fn to_ad8(self) -> AD8 {
//         self.into()
//     }
//
//     fn to_ad9(self) -> AD9 {
//         self.into()
//     }
//
//     fn to_ad10(self) -> AD10 {
//         self.into()
//     }
// }

//impl AD for f64 {}

/// Lift AD functions
///
/// # Description
/// To lift `AD` functions
///
/// # Implementation
///
/// * All `Fn(T) -> T where T:AD` functions can be lift to `Fn(f64) -> f64`
/// * If `j > i`, then `Fn(AD{j}) -> AD{j}` can be lift to `Fn(AD{i}) -> AD{i}`
///
/// # Usage
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let ad0 = 2f64;
///     let ad1 = AD1::new(2f64, 1f64);
///     
///     let lift1 = ADLift::<_, AD1>::new(f_ad);
///     let lift2 = ADLift::<_, AD2>::new(f_ad2);
///
///     let ans_ad0 = ad0.powi(2);
///     let ans_ad1 = ad1.powi(2);
///
///     // All AD function can be lift to f64
///     assert_eq!(ans_ad0, lift1.call_stable(ad0));
///     assert_eq!(ans_ad0, lift2.call_stable(ad0));
///
///     // AD2 is higher than AD1 (AD2 -> AD1 lifting is allowed)
///     assert_eq!(ans_ad1, lift2.call_stable(ad1));
/// }
///
/// fn f_ad<T: AD>(x: T) -> T {
///     x.powi(2)
/// }
///
/// fn f_ad2(x: AD2) -> AD2 {
///     x.powi(2)
/// }
/// ```
pub struct ADLift<F, T> {
    f: Box<F>,
    _marker: PhantomData<T>,
}

impl<F: Fn(T) -> T, T> ADLift<F, T> {
    pub fn new(f: F) -> Self {
        Self {
            f: Box::new(f),
            _marker: PhantomData
        }
    }

    pub fn f(&self, t: T) -> T {
        (self.f)(t)
    }
}

impl<F: Fn(T) -> T, T: AD> StableFn<f64> for ADLift<F, T> {
    type Output = f64;
    fn call_stable(&self, target: f64) -> Self::Output {
        self.f(T::from(target)).into()
    }
}

impl<F: Fn(T) -> T, T: AD> StableFn<AD1> for ADLift<F, T> {
    type Output = AD1;
    fn call_stable(&self, target: AD1) -> Self::Output {
        self.f(T::from(target)).into()
    }
}

// Nightly only
//pub struct ADLift<F>(F);
//
//impl<F: FnOnce<(T)>, T> FnOnce<f64> for ADLift<F> {
//    type Output = f64;
//    
//}

