//! Dual number system for Automatic differentiation
//!
//! ## Dual number system
//!
//! * `Dual` is structure for AD
//!     * `value(&self) -> f64` : Value
//!     * `slope(&self) -> f64` : Slope (value of derivatives)
//!     * `extract(&self) -> (f64,f64)` : Extract both
//! * Constructor for `Dual`
//!     * `Dual::new(T, T)`
//!     * `dual(T, T)`
//! * Implemented Operations (`Real` trait is implemented)
//!     * `Add, Sub, Mul, Div`
//!     * `sin, cos, tan`
//!     * `asin, acos, atan`
//!     * `sinh, cosh, tanh`
//!     * `asinh, acosh, atanh`
//!     * `sin_cos`
//!     * `exp, ln, log, log2, log10`
//!     * `powi, powf, sqrt`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let x = dual(1, 1); // x at x = 1
//!         (x.clone() + x.clone()).print(); // dual(2, 2)
//!         (x.clone() - x.clone()).print(); // dual(0, 0)
//!         (x.clone() * x.clone()).print(); // dual(1, 2)
//!         (x.clone() / x.clone()).print(); // dual(1, 0)
//!         x.sin().print();
//!         x.cos().print();
//!         x.exp().print();
//!         x.ln().print();
//!         x.powi(2).print(); // same as x * x
//!     }
//!     ```
//!
//!     * After ver `0.10.1`, you can use reference operations
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let x = dual(1,1);
//!         (&x + &x).print();
//!         (&x - &x).print();
//!         // and etc.
//!     }
//!     ```

use operation::extra_ops::{ExpLogOps, PowOps, TrigOps};
use std::convert::Into;
/// Structure for Automatic Differentiation
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};
use structure::vector::*;
//use self::num_traits::{Num, Zero, One, NumCast, ToPrimitive};
//use std::num::ParseFloatError;

/// Dual - Structure for AD
///
/// # Fields
/// `x: f64` : value
/// `dx: f64` : slope at `x`
///
/// But they are private fields.
/// You should use `value` or `slope` function to extract them.
#[derive(Debug, Copy, Clone, Default)]
pub struct Dual {
    x: f64,
    dx: f64,
}

impl fmt::Display for Dual {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s1 = format!("value: {}\nslope: {}", self.x, self.dx);
        write!(f, "{}", s1)
    }
}

impl PartialEq for Dual {
    fn eq(&self, other: &Dual) -> bool {
        self.x == other.x && self.dx == other.dx
    }
}

impl Dual {
    /// Syntactic sugar to declare dual
    ///
    /// # Type
    /// `(T, T) -> Dual where T: Into<f64> + Copy`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = Dual::new(2, 1); // y = x at x = 2
    /// a.print();
    /// ```
    pub fn new<T: Into<f64> + Copy>(x: T, dx: T) -> Dual {
        Dual {
            x: x.into(),
            dx: dx.into(),
        }
    }

    /// Just extract value
    pub fn value(&self) -> f64 {
        self.x
    }

    /// Just extract slope
    pub fn slope(&self) -> f64 {
        self.dx
    }

    /// Just extract both
    pub fn extract(&self) -> (f64, f64) {
        (self.x, self.dx)
    }
}

pub fn dual<T: Into<f64> + Copy>(x: T, dx: T) -> Dual {
    Dual::new(x, dx)
}

// =============================================================================
// Ops with each other
// =============================================================================

/// Neg for Dual
impl Neg for Dual {
    type Output = Dual;
    fn neg(self) -> Dual {
        Dual::new(-self.x, -self.dx)
    }
}

/// Add for Dual
impl Add<Dual> for Dual {
    type Output = Dual;
    fn add(self, other: Dual) -> Dual {
        Dual::new(self.x + other.x, self.dx + other.dx)
    }
}

/// Sub for Dual
impl Sub<Dual> for Dual {
    type Output = Dual;
    fn sub(self, other: Dual) -> Dual {
        Dual::new(self.x - other.x, self.dx - other.dx)
    }
}

/// Mul for Dual
impl Mul<Dual> for Dual {
    type Output = Dual;
    fn mul(self, other: Dual) -> Dual {
        let v1 = self.x;
        let v2 = other.x;
        let dv1 = self.dx;
        let dv2 = other.dx;

        Dual::new(v1 * v2, v1 * dv2 + v2 * dv1)
    }
}

/// Div for Dual
impl Div<Dual> for Dual {
    type Output = Dual;
    fn div(self, other: Dual) -> Dual {
        let v1 = self.x;
        let v2 = other.x;
        let dv1 = self.dx;
        let dv2 = other.dx;

        Dual::new(v1 / v2, (dv1 * v2 - v1 * dv2) / (v2 * v2))
    }
}

impl Rem<Dual> for Dual {
    type Output = Dual;
    fn rem(self, _rhs: Dual) -> Dual {
        unimplemented!()
    }
}

impl<'a> Neg for &'a Dual {
    type Output = Dual;

    fn neg(self) -> Self::Output {
        Dual {
            x: -self.x,
            dx: -self.dx,
        }
    }
}

impl<'a, 'b> Add<&'b Dual> for &'a Dual {
    type Output = Dual;

    fn add(self, rhs: &Dual) -> Self::Output {
        Dual {
            x: self.x + rhs.x,
            dx: self.dx + rhs.dx,
        }
    }
}

impl<'a, 'b> Sub<&'b Dual> for &'a Dual {
    type Output = Dual;

    fn sub(self, rhs: &Dual) -> Self::Output {
        Dual {
            x: self.x - rhs.x,
            dx: self.dx - rhs.dx,
        }
    }
}

impl<'a, 'b> Mul<&'b Dual> for &'a Dual {
    type Output = Dual;

    fn mul(self, rhs: &Dual) -> Self::Output {
        Dual {
            x: self.x * rhs.x,
            dx: self.x * rhs.dx + self.dx * rhs.x,
        }
    }
}

impl<'a, 'b> Div<&'b Dual> for &'a Dual {
    type Output = Dual;

    fn div(self, rhs: &Dual) -> Self::Output {
        Dual {
            x: self.x / rhs.x,
            dx: (self.dx * rhs.x - self.x * rhs.dx) / (rhs.x * rhs.x),
        }
    }
}

// =============================================================================
// Ops with f64
// =============================================================================

impl Add<f64> for Dual {
    type Output = Dual;
    fn add(self, other: f64) -> Dual {
        self + Dual::new(other, 0.)
    }
}

impl Sub<f64> for Dual {
    type Output = Dual;
    fn sub(self, other: f64) -> Dual {
        self - Dual::new(other, 0.)
    }
}

impl Mul<f64> for Dual {
    type Output = Dual;
    fn mul(self, other: f64) -> Dual {
        self * Dual::new(other, 0.)
    }
}

impl Div<f64> for Dual {
    type Output = Dual;
    fn div(self, other: f64) -> Dual {
        self / Dual::new(other, 0.)
    }
}

impl Add<Dual> for f64 {
    type Output = Dual;
    fn add(self, other: Dual) -> Dual {
        other.add(self)
    }
}

impl Sub<Dual> for f64 {
    type Output = Dual;
    fn sub(self, other: Dual) -> Dual {
        dual(self, 0.) - other
    }
}

impl Mul<Dual> for f64 {
    type Output = Dual;
    fn mul(self, other: Dual) -> Dual {
        other.mul(self)
    }
}

impl Div<Dual> for f64 {
    type Output = Dual;
    fn div(self, other: Dual) -> Dual {
        dual(self, 0.) / other
    }
}

impl<'a> Add<f64> for &'a Dual {
    type Output = Dual;

    fn add(self, rhs: f64) -> Self::Output {
        Dual {
            x: self.x + rhs,
            dx: self.dx,
        }
    }
}

impl<'a> Sub<f64> for &'a Dual {
    type Output = Dual;

    fn sub(self, rhs: f64) -> Self::Output {
        Dual {
            x: self.x - rhs,
            dx: self.dx,
        }
    }
}

impl<'a> Mul<f64> for &'a Dual {
    type Output = Dual;

    fn mul(self, rhs: f64) -> Self::Output {
        Dual {
            x: self.x * rhs,
            dx: self.dx * rhs,
        }
    }
}

impl<'a> Div<f64> for &'a Dual {
    type Output = Dual;

    fn div(self, rhs: f64) -> Self::Output {
        Dual {
            x: self.x / rhs,
            dx: self.dx / rhs,
        }
    }
}

// =============================================================================
// Trigonometric Ops
// =============================================================================

/// Trigonometric function with Dual
impl TrigOps for Dual {
    fn sin(&self) -> Self {
        let val = self.x.sin();
        let dval = self.dx * self.x.cos();
        Dual::new(val, dval)
    }

    fn cos(&self) -> Self {
        let val = self.x.cos();
        let dval = -self.dx * self.x.sin();
        Dual::new(val, dval)
    }

    fn tan(&self) -> Self {
        let val = self.x.tan();
        let dval = self.dx * (1. + val * val); // 1 + tan^2 = sec^2
        Dual::new(val, dval)
    }

    fn asin(&self) -> Self {
        let val = self.x.asin();
        let dval = self.dx / (1f64 - self.x.powi(2)).sqrt();
        Dual::new(val, dval)
    }

    fn acos(&self) -> Self {
        let val = self.x.acos();
        let dval = -self.dx / (1f64 - self.x.powi(2)).sqrt();
        Dual::new(val, dval)
    }

    fn atan(&self) -> Self {
        let val = self.x.atan();
        let dval = self.dx / (1f64 + self.x.powi(2));
        Dual::new(val, dval)
    }

    fn sinh(&self) -> Self {
        let val = self.x.sinh();
        let dval = self.dx * self.x.cosh();
        Dual::new(val, dval)
    }

    fn cosh(&self) -> Self {
        let val = self.x.cosh();
        let dval = self.dx * self.x.sinh();
        Dual::new(val, dval)
    }

    fn tanh(&self) -> Self {
        let val = self.x.tanh();
        let dval = self.dx / self.x.cosh().powi(2);
        Dual::new(val, dval)
    }

    fn asinh(&self) -> Self {
        let val = self.x.asinh();
        let dval = self.dx / (1f64 + self.x.powi(2)).sqrt();
        Dual::new(val, dval)
    }

    fn acosh(&self) -> Self {
        let val = self.x.acosh();
        let dval = self.dx / (self.x.powi(2) - 1f64).sqrt();
        Dual::new(val, dval)
    }

    fn atanh(&self) -> Self {
        let val = self.x.atanh();
        let dval = self.dx / (1f64 - self.x.powi(2));
        Dual::new(val, dval)
    }

    fn sin_cos(&self) -> (Self, Self) {
        let vals = self.x.sin_cos();
        let dvals = (self.dx * vals.1, -self.dx * vals.0);
        let d1 = Dual::new(vals.0, dvals.0);
        let d2 = Dual::new(vals.1, dvals.1);
        (d1, d2)
    }
}

// =============================================================================
// Exp & Logarithm
// =============================================================================

/// Exponential & Logarithm for Dual
impl ExpLogOps for Dual {
    fn exp(&self) -> Self {
        let val = self.value().exp();
        let dval = val * self.slope();
        Dual::new(val, dval)
    }

    fn ln(&self) -> Self {
        assert_ne!(self.value(), 0.);
        let val = self.value().ln();
        let dval = self.slope() / self.value();
        Dual::new(val, dval)
    }

    fn log(&self, base: f64) -> Self {
        self.ln() / base.ln()
    }

    fn log2(&self) -> Self {
        self.log(2f64)
    }

    fn log10(&self) -> Self {
        self.log(10f64)
    }
}

// =============================================================================
// Power
// =============================================================================
/// Power for Dual
impl PowOps for Dual {
    fn powi(&self, n: i32) -> Self {
        let x = self.x;
        let val = x.powi(n);
        let dval = (n as f64) * x.powi(n - 1) * self.dx;
        Dual::new(val, dval)
    }

    fn powf(&self, f: Self) -> Self {
        let x = self.value();
        let dx = self.slope();
        let n = f.value();
        let dn = f.slope();
        Dual::new(x.powf(n), x.powf(n - 1f64) * (n * dx + x * x.ln() * dn))
    }

    fn sqrt(&self) -> Self {
        self.powf(Dual::new(0.5f64, 0f64))
    }
}

// =============================================================================
// Dual List
// =============================================================================
/// Convert Vector <=> Dual
pub trait VecWithDual {
    type Item;
    fn conv_dual(&self) -> Self::Item;
}

impl VecWithDual for Vec<f64> {
    type Item = Vec<Dual>;
    fn conv_dual(&self) -> Vec<Dual> {
        self.clone()
            .into_iter()
            .map(|x| Dual::new(x, 0.))
            .collect::<Vec<Dual>>()
    }
}

impl VecWithDual for Vec<Dual> {
    type Item = Vec<f64>;
    fn conv_dual(&self) -> Vec<f64> {
        self.clone()
            .into_iter()
            .map(|x| x.value())
            .collect::<Vec<f64>>()
    }
}

pub trait Dualist {
    fn values(&self) -> Vec<f64>;
    fn slopes(&self) -> Vec<f64>;
}

impl Dualist for Vec<Dual> {
    fn values(&self) -> Vec<f64> {
        self.conv_dual()
    }

    fn slopes(&self) -> Vec<f64> {
        self.clone()
            .into_iter()
            .map(|t| t.slope())
            .collect::<Vec<f64>>()
    }
}

pub fn merge_dual(y: Vec<f64>, dy: Vec<f64>) -> Vec<Dual> {
    y.into_iter()
        .zip(&dy)
        .map(|(t, &dt)| Dual::new(t, dt))
        .collect::<Vec<Dual>>()
}

impl FPVector for Vec<Dual> {
    type Scalar = Dual;

    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> Self::Scalar,
    {
        self.clone()
            .into_iter()
            .map(|x| f(x))
            .collect::<Vec<Dual>>()
    }

    fn reduce<F, T>(&self, _init: T, _f: F) -> Self::Scalar
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
        T: Into<Self::Scalar>,
    {
        unimplemented!()
    }

    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
    {
        self.into_iter()
            .zip(other)
            .map(|(x, y)| f(*x, *y))
            .collect::<Vec<Dual>>()
    }

    fn filter<F>(&self, _f: F) -> Self
    where
        F: Fn(Self::Scalar) -> bool,
    {
        unimplemented!()
    }

    fn take(&self, n: usize) -> Self {
        self.clone().into_iter().take(n).collect::<Vec<Dual>>()
    }

    fn skip(&self, n: usize) -> Self {
        self.clone().into_iter().skip(n).collect::<Vec<Dual>>()
    }
}

#[allow(unused_variables)]
impl VecOps for Vec<Dual> {
    type Scalar = Dual;

    fn add(&self, other: &Self) -> Self {
        self.zip_with(|x, y| x + y, other)
    }

    fn sub(&self, other: &Self) -> Self {
        self.zip_with(|x, y| x - y, other)
    }

    fn mul(&self, other: &Self) -> Self {
        self.zip_with(|x, y| x * y, other)
    }

    fn div(&self, other: &Self) -> Self {
        self.zip_with(|x, y| x / y, other)
    }

    fn dot(&self, _other: &Self) -> Self::Scalar {
        // dot product of Dual is similar to Complex with \epsilon^2 = 0
        self.clone()
            .into_iter()
            .zip(_other.clone())
            .map(|(x, y)| x.mul(y))
            .fold(Self::Scalar::new(0., 0.), |sum: Self::Scalar, x| sum.add(x))
    }

    fn norm(&self) -> Self::Scalar {
        unimplemented!()
    }

    fn normalize(&self) -> Self {
        unimplemented!()
    }

    fn s_add(&self, scala: f64) -> Self {
        self.fmap(|x| x + scala)
    }

    fn s_sub(&self, scala: f64) -> Self {
        self.fmap(|x| x - scala)
    }

    fn s_mul(&self, scala: f64) -> Self {
        self.fmap(|x| x * scala)
    }

    fn s_div(&self, scala: f64) -> Self {
        self.fmap(|x| x / scala)
    }
}

// =============================================================================
// Num Traits for Dual
// =============================================================================
//impl Zero for Dual {
//    fn zero() -> Self {
//        Dual::new(0f64, 0f64)
//    }
//
//    fn is_zero(&self) -> bool {
//        self.x == 0f64 && self.dx == 0f64
//    }
//}
//
//impl One for Dual {
//    fn one() -> Self {
//        Dual::new(1f64, 0f64)
//    }
//}
//
//impl Num for Dual {
//    type FromStrRadixErr = ParseFloatError;
//
//    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
//        unimplemented!()
//    }
//}
//
//impl ToPrimitive for Dual {
//    fn to_i64(&self) -> Option<i64> {
//        Some(self.x as i64)
//    }
//
//    fn to_u64(&self) -> Option<u64> {
//        Some(self.x as u64)
//    }
//
//    fn to_f64(&self) -> Option<f64> {
//        Some(self.x)
//    }
//}
//
//impl NumCast for Dual {
//    fn from<T: ToPrimitive>(n: T) -> Option<Self> {
//        match n.to_f64() {
//            None => None,
//            Some(val) => Some(Dual::new(val, 0f64))
//        }
//    }
//}
//
//impl Real for Dual {
//    fn min_value() -> Self {
//        unimplemented!()
//    }
//
//    fn min_positive_value() -> Self {
//        unimplemented!()
//    }
//
//    fn epsilon() -> Self {
//        unimplemented!()
//    }
//
//    fn max_value() -> Self {
//        unimplemented!()
//    }
//
//    fn floor(self) -> Self {
//        unimplemented!()
//    }
//
//    fn ceil(self) -> Self {
//        unimplemented!()
//    }
//
//    fn round(self) -> Self {
//        unimplemented!()
//    }
//
//    fn trunc(self) -> Self {
//        unimplemented!()
//    }
//
//    fn fract(self) -> Self {
//        unimplemented!()
//    }
//
//    fn abs(self) -> Self {
//        unimplemented!()
//    }
//
//    fn signum(self) -> Self {
//        unimplemented!()
//    }
//
//    fn is_sign_positive(self) -> bool {
//        unimplemented!()
//    }
//
//    fn is_sign_negative(self) -> bool {
//        unimplemented!()
//    }
//
//    fn mul_add(self, a: Self, b: Self) -> Self {
//        unimplemented!()
//    }
//
//    fn recip(self) -> Self {
//        let val = 1f64 / self.x;
//        let dval = - self.dx / self.x.powi(2);
//    }
//
//    fn powi(self, n: i32) -> Self {
//        let x = self.x;
//        let val = x.powi(n);
//        let dval = (n as f64) * x.powi(n - 1) * self.dx;
//        Dual::new(val, dval)
//    }
//
//    fn powf(self, n: Self) -> Self {
//        let f = n.to_f64();
//        let val = self.x.powf(f);
//        let dval = f * self.x.powf(f - 1f64) * self.dx;
//        Dual::new(val, dval)
//    }
//
//    fn sqrt(self) -> Self {
//        self.powf(Dual::new(0.5f64, 0f64))
//    }
//
//    fn exp(self) -> Self {
//        let val = self.x.exp();
//        let dval = val * self.dx;
//        Dual::new(val, dval)
//    }
//
//    fn exp2(self) -> Self {
//        let val = self.x.exp2();
//        let dval = val * 2f64.ln() * self.dx;
//        Dual::new(val, dval)
//    }
//
//    fn ln(self) -> Self {
//        let val = self.x.ln();
//        let dval = self.dx / self.x;
//        Dual::new(val, dval)
//    }
//
//    fn log(self, base: Self) -> Self {
//        self.ln() / base.to_f64().ln()
//    }
//
//    fn log2(self) -> Self {
//        self.log(Dual::new(2f64, 0f64))
//    }
//
//    fn log10(self) -> Self {
//        self.log(Dual::new(10f64, 0f64))
//    }
//
//    fn to_degrees(self) -> Self {
//        unimplemented!()
//    }
//
//    fn to_radians(self) -> Self {
//        unimplemented!()
//    }
//
//    fn max(self, other: Self) -> Self {
//        unimplemented!()
//    }
//
//    fn min(self, other: Self) -> Self {
//        unimplemented!()
//    }
//
//    fn abs_sub(self, other: Self) -> Self {
//        unimplemented!()
//    }
//
//    fn cbrt(self) -> Self {
//        unimplemented!()
//    }
//
//    fn hypot(self, other: Self) -> Self {
//        unimplemented!()
//    }
//
//    fn sin(self) -> Self {
//        let val = self.x.sin();
//        let dval = self.dx * self.x.cos();
//        Dual::new(val, dval)
//    }
//
//    fn cos(self) -> Self {
//        let val = self.x.cos();
//        let dval = -self.dx * self.x.sin();
//        Dual::new(val, dval)
//    }
//
//    fn tan(self) -> Self {
//        let val = self.x.tan();
//        let dval = self.dx * (1. + val * val); // 1 + tan^2 = sec^2
//        Dual::new(val, dval)
//    }
//
//    fn asin(self) -> Self {
//        let val = self.x.asin();
//        let dval = self.dx / (1f64 - self.x.powi(2)).sqrt();
//        Dual::new(val, dval)
//    }
//
//    fn acos(self) -> Self {
//        let val = self.x.acos();
//        let dval = - self.dx / (1f64 - self.x.powi(2)).sqrt();
//        Dual::new(val, dval)
//    }
//
//    fn atan(self) -> Self {
//        let val = self.x.atan();
//        let dval = self.dx / (1f64 + self.x.powi(2));
//        Dual::new(val, dval)
//    }
//
//    fn atan2(self, other: Self) -> Self {
//        unimplemented!()
//    }
//
//    fn sin_cos(self) -> (Self, Self) {
//        let vals = self.sin_cos();
//        let dvals = (self.dx * vals.1, -self.dx * vals.0);
//        (Dual::new(vals.0, dvals.0), Dual::new(vals.1, dvals.1))
//    }
//
//    fn exp_m1(self) -> Self {
//        self.exp() - 1f64
//    }
//
//    fn ln_1p(self) -> Self {
//        unimplemented!()
//    }
//
//    fn sinh(self) -> Self {
//        let val = self.x.sinh();
//        let dval = self.dx * self.x.cosh();
//        Dual::new(val, dval)
//    }
//
//    fn cosh(self) -> Self {
//        let val = self.x.cosh();
//        let dval = self.dx * self.x.sinh();
//        Dual::new(val, dval)
//    }
//
//    fn tanh(self) -> Self {
//        let val = self.x.tanh();
//        let dval = self.dx / self.x.cosh().powi(2);
//        Dual::new(val, dval)
//
//    }
//
//    fn asinh(self) -> Self {
//        let val = self.x.asinh();
//        let dval = self.dx / (1f64 + self.x.powi(2)).sqrt();
//        Dual::new(val, dval)
//    }
//
//    fn acosh(self) -> Self {
//        let val = self.x.acosh();
//        let dval = self.dx / (self.x.powi(2) - 1f64).sqrt();
//        Dual::new(val, dval)
//    }
//
//    fn atanh(self) -> Self {
//        let val = self.x.atanh();
//        let dval = self.dx / (1f64 - self.x.powi(2));
//        Dual::new(val, dval)
//    }
//}
