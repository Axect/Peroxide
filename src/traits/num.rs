//! Missing operations & comprehensive number structures
//!
//! ## `Real` trait
//!
//! * `Real` is a trait for binding `f64`, `Dual`, `HyperDual`
//! * `Real` requires `PowOps, TrigOps, ExpLogOps` & `std::Ops<Self>` & `std::Ops<f64>`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
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
use crate::structure::{
    dual::Dual,
    hyper_dual::HyperDual,
};
use self::Number::{F, D};

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
// Number
// =============================================================================

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub enum Number {
    F(f64),
    D(Dual),
}

impl Default for Number {
    fn default() -> Self {
        F(f64::default())
    }
}

impl Neg for Number {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            F(x) => F(-x),
            D(x) => D(-x),
        }
    }
}

impl Add for Number {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (F(x), F(y)) => F(x + y),
            (D(x), D(y)) => D(x + y),
            (F(x), D(y)) => D(x + y),
            (D(x), F(y)) => D(x + y),
        }
    }
}

impl Add<f64> for Number {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        self.add(F(rhs))
    }
}

impl Add<Dual> for Number {
    type Output = Self;

    fn add(self, rhs: Dual) -> Self::Output {
        self.add(D(rhs))
    }
}

impl Sub for Number {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (F(x), F(y)) => F(x - y),
            (D(x), D(y)) => D(x - y),
            (F(x), D(y)) => D(x - y),
            (D(x), F(y)) => D(x - y),
        }
    }
}

impl Sub<f64> for Number {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        self.sub(F(rhs))
    }
}

impl Sub<Dual> for Number {
    type Output = Self;

    fn sub(self, rhs: Dual) -> Self::Output {
        self.sub(D(rhs))
    }
}

impl Mul for Number {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (F(x), F(y)) => F(x * y),
            (D(x), D(y)) => D(x * y),
            (F(x), D(y)) => D(x * y),
            (D(x), F(y)) => D(x * y),
        }
    }
}

impl Mul<f64> for Number {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        self.mul(F(rhs))
    }
}

impl Mul<Dual> for Number {
    type Output = Self;

    fn mul(self, rhs: Dual) -> Self::Output {
        self.mul(D(rhs))
    }
}

impl Div for Number {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (F(x), F(y)) => F(x / y),
            (D(x), D(y)) => D(x / y),
            (F(x), D(y)) => D(x / y),
            (D(x), F(y)) => D(x / y),
        }
    }
}

impl Div<f64> for Number {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        self.div(F(rhs))
    }
}

impl Div<Dual> for Number {
    type Output = Self;

    fn div(self, rhs: Dual) -> Self::Output {
        self.div(D(rhs))
    }
}

impl ExpLogOps for Number {
    fn exp(&self) -> Self {
        match self {
            F(x) => F(x.exp()),
            D(x) => D(x.exp()),
        }
    }

    fn ln(&self) -> Self {
        match self {
            F(x) => F(x.ln()),
            D(x) => D(x.exp()),
        }
    }

    fn log(&self, base: f64) -> Self {
        match self {
            F(x) => F(x.log(base)),
            D(x) => D(x.log(base)),
        }
    }

    fn log2(&self) -> Self {
        match self {
            F(x) => F(x.log2()),
            D(x) => D(x.log2()),
        }
    }

    fn log10(&self) -> Self {
        match self {
            F(x) => F(x.log10()),
            D(x) => D(x.log10()),
        }
    }
}

impl TrigOps for Number {
    fn sin(&self) -> Self {
        match self {
            F(x) => F(x.sin()),
            D(x) => D(x.sin()),
        }
    }

    fn cos(&self) -> Self {
        match self {
            F(x) => F(x.cos()),
            D(x) => D(x.cos()),
        }
    }

    fn tan(&self) -> Self {
        match self {
            F(x) => F(x.tan()),
            D(x) => D(x.tan()),
        }
    }

    fn asin(&self) -> Self {
        match self {
            F(x) => F(x.asin()),
            D(x) => D(x.asin()),
        }
    }

    fn acos(&self) -> Self {
        match self {
            F(x) => F(x.acos()),
            D(x) => D(x.acos()),
        }
    }

    fn atan(&self) -> Self {
        match self {
            F(x) => F(x.atan()),
            D(x) => D(x.atan()),
        }
    }

    fn sinh(&self) -> Self {
        match self {
            F(x) => F(x.sinh()),
            D(x) => D(x.sinh()),
        }
    }

    fn cosh(&self) -> Self {
        match self {
            F(x) => F(x.cosh()),
            D(x) => D(x.cosh()),
        }
    }

    fn tanh(&self) -> Self {
        match self {
            F(x) => F(x.tanh()),
            D(x) => D(x.tanh()),
        }
    }

    fn asinh(&self) -> Self {
        match self {
            F(x) => F(x.asinh()),
            D(x) => D(x.asinh()),
        }
    }

    fn acosh(&self) -> Self {
        match self {
            F(x) => F(x.acosh()),
            D(x) => D(x.acosh()),
        }
    }

    fn atanh(&self) -> Self {
        match self {
            F(x) => F(x.atanh()),
            D(x) => D(x.atanh()),
        }
    }

    fn sin_cos(&self) -> (Self, Self) {
        match self {
            F(x) => (F(x.sin()), F(x.cos())),
            D(x) => (D(x.sin()), D(x.cos())),
        }
    }
}

impl PowOps for Number {
    fn powi(&self, n: i32) -> Self {
        match self {
            F(x) => F(x.powi(n)),
            D(x) => D(x.powi(n)),
        }
    }

    fn powf(&self, f: f64) -> Self {
        match self {
            F(x) => F(x.powf(f)),
            D(x) => D(x.powf(f)),
        }
    }

    #[allow(unreachable_patterns)]
    fn pow(&self, f: Self) -> Self {
        match (self, &f) {
            (F(x), F(y)) => F(x.powf(*y)),
            (F(x), D(y)) => D(Dual::from_f64(*x).pow(*y)),
            (D(x), F(y)) => D(x.pow(Dual::from_f64(*y))),
            (D(x), D(y)) => D(x.pow(*y)),
            _ => unreachable!(),
        }
    }

    fn sqrt(&self) -> Self {
        match self {
            F(x) => F(x.sqrt()),
            D(x) => D(x.sqrt()),
        }
    }
}

impl Real for Number {
    fn to_f64(&self) -> f64 {
        match self {
            F(x) => x.to_owned(),
            D(x) => x.to_f64(),
        }
    }

    fn from_f64(f: f64) -> Self {
        F(f)
    }

    fn to_dual(&self) -> Dual {
        match self {
            F(x) => x.to_dual(),
            D(x) => x.to_owned(),
        }
    }

    fn from_dual(d: Dual) -> Self {
        D(d)
    }

    fn to_hyper_dual(&self) -> HyperDual {
        unimplemented!()
    }

    fn from_hyper_dual(_h: HyperDual) -> Self {
        unimplemented!()
    }
}

impl Add<Number> for f64 {
    type Output = Number;

    fn add(self, rhs: Number) -> Self::Output {
        rhs.add(self)
    }
}

impl Sub<Number> for f64 {
    type Output = Number;

    fn sub(self, rhs: Number) -> Self::Output {
        -rhs.sub(self)
    }
}

impl Mul<Number> for f64 {
    type Output = Number;

    fn mul(self, rhs: Number) -> Self::Output {
        rhs.mul(self)
    }
}

impl Div<Number> for f64 {
    type Output = Number;

    fn div(self, rhs: Number) -> Self::Output {
        match rhs {
            F(x) => F(self / x),
            D(x) => D(self / x),
        }
    }
}

pub trait NumberVector {
    fn to_dual_vec(&self) -> Vec<Dual>;
    fn to_f64_vec(&self) -> Vec<f64>;
    fn to_hyper_vec(&self) -> Vec<HyperDual>;
    fn from_dual_vec(v: Vec<Dual>) -> Self;
    fn from_f64_vec(v: Vec<f64>) -> Self;
    fn from_hyper_vec(v: Vec<HyperDual>) -> Self;
}

impl NumberVector for Vec<Number> {
    fn to_dual_vec(&self) -> Vec<Dual> {
        self.clone().into_iter().map(|x| x.to_dual()).collect()
    }

    fn to_f64_vec(&self) -> Vec<f64> {
        self.clone().into_iter().map(|x| x.to_f64()).collect()
    }

    fn to_hyper_vec(&self) -> Vec<HyperDual> {
        self.clone()
            .into_iter()
            .map(|x| x.to_hyper_dual())
            .collect()
    }

    fn from_dual_vec(v: Vec<Dual>) -> Self {
        v.into_iter().map(|x| Number::from_dual(x)).collect()
    }

    fn from_f64_vec(v: Vec<f64>) -> Self {
        v.into_iter().map(|x| Number::from_f64(x)).collect()
    }

    fn from_hyper_vec(v: Vec<HyperDual>) -> Self {
        v.into_iter().map(|x| Number::from_hyper_dual(x)).collect()
    }
}
