//! Missing operations & comprehensive number structures
//!
//! ## `Numeric` trait
//!
//! * `Numeric` requires `PowOps, TrigOps, ExpLogOps` & `std::Ops<Self>` & `std::Ops<f64>`

use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait PowOps: Sized {
    type Float;
    fn powi(&self, n: i32) -> Self;
    fn powf(&self, f: Self::Float) -> Self;
    fn pow(&self, f: Self) -> Self;
    fn sqrt(&self) -> Self;
}

pub trait TrigOps: Sized + Div<Output = Self> {
    fn sin_cos(&self) -> (Self, Self);
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
    fn sinh(&self) -> Self;
    fn cosh(&self) -> Self;
    fn tanh(&self) -> Self;
    fn asin(&self) -> Self;
    fn acos(&self) -> Self;
    fn atan(&self) -> Self;
    fn asin_acos(&self) -> (Self, Self) {
        (self.asin(), self.acos())
    }
    fn asinh(&self) -> Self;
    fn acosh(&self) -> Self;
    fn atanh(&self) -> Self;
    fn asinh_acosh(&self) -> (Self, Self) {
        (self.asinh(), self.acosh())
    }
}

pub trait ExpLogOps: Sized {
    type Float;
    fn exp(&self) -> Self;
    fn ln(&self) -> Self;
    fn log(&self, base: Self::Float) -> Self;
    fn log2(&self) -> Self;
    fn log10(&self) -> Self;
}

pub trait Float:
    PowOps<Float = Self>
    + TrigOps
    + ExpLogOps<Float = Self>
    + Neg<Output = Self>
    + Add<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Sub<Output = Self>
    + PartialOrd
    + Copy
    + Clone
{
}

pub trait Numeric<T: Float>:
    PowOps<Float = T>
    + TrigOps
    + ExpLogOps<Float = T>
    + Neg<Output = Self>
    + Add<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Sub<Output = Self>
    + Add<T, Output = Self>
    + Mul<T, Output = Self>
    + Div<T, Output = Self>
    + Sub<T, Output = Self>
    + Clone
{
}

// =============================================================================
// Numeric Traits for f32, f64
// =============================================================================

macro_rules! impl_float {
    ($($t:ty),*) => {
        $(
            impl PowOps for $t {
                type Float = $t;
                fn powi(&self, n: i32) -> Self {
                    (*self).powi(n)
                }

                fn powf(&self, f: Self::Float) -> Self {
                    (*self).powf(f)
                }

                fn pow(&self, f: Self) -> Self {
                    (*self).powf(f)
                }

                fn sqrt(&self) -> Self {
                    (*self).sqrt()
                }
            }

            impl TrigOps for $t {
                fn sin_cos(&self) -> (Self, Self) {
                    (*self).sin_cos()
                }

                fn sin(&self) -> Self {
                    (*self).sin()
                }

                fn cos(&self) -> Self {
                    (*self).cos()
                }

                fn tan(&self) -> Self {
                    (*self).tan()
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

                fn asin(&self) -> Self {
                    (*self).asin()
                }

                fn acos(&self) -> Self {
                    (*self).acos()
                }

                fn atan(&self) -> Self {
                    (*self).atan()
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
            }

            impl ExpLogOps for $t {
                type Float = $t;
                fn exp(&self) -> Self {
                    (*self).exp()
                }
                fn ln(&self) -> Self {
                    (*self).ln()
                }
                fn log(&self, base: Self::Float) -> Self {
                    (*self).log(base)
                }
                fn log2(&self) -> Self {
                    (*self).log2()
                }
                fn log10(&self) -> Self {
                    (*self).log10()
                }
            }

            impl Float for $t {}
        )*
    };
}

impl_float!(f32, f64);

impl Numeric<f32> for f32 {}
impl Numeric<f64> for f64 {}

// ┌──────────────────────────────────────────────────────────┐
//  Fundamental Traits
// └──────────────────────────────────────────────────────────┘
pub trait Group: Sized + Add<Self, Output = Self> {
    fn zero() -> Self;
}

pub trait Ring: Group + Mul<Self, Output = Self> {
    fn one() -> Self;
}

impl Group for f32 {
    fn zero() -> Self {
        0.0
    }
}

impl Ring for f32 {
    fn one() -> Self {
        1.0
    }
}

impl Group for f64 {
    fn zero() -> Self {
        0.0
    }
}

impl Ring for f64 {
    fn one() -> Self {
        1.0
    }
}
