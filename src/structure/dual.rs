use operation::extra_ops::{ExpLogOps, PowOps, TrigOps};
use std::convert::Into;
/// Structure for Automatic Differentiation
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};
use structure::vector::*;
use std::f64::consts::E;

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
        assert_ne!(other.x, 0f64);
        let v1 = self.x;
        let v2 = other.x;
        let dv1 = self.dx;
        let dv2 = other.dx;

        Dual::new(v1 / v2, (dv1 * v2 - v1 * dv2) / (v2 * v2))
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

// =============================================================================
// Trigonometric Ops
// =============================================================================

/// Trigonometric function with Dual
impl TrigOps for Dual {
    type Output = Dual;

    fn sin(&self) -> Dual {
        let val = self.x.sin();
        let dval = self.dx * self.x.cos();
        Dual::new(val, dval)
    }

    fn cos(&self) -> Dual {
        let val = self.x.cos();
        let dval = -self.dx * self.x.sin();
        Dual::new(val, dval)
    }

    fn tan(&self) -> Dual {
        let val = self.x.tan();
        let dval = self.dx * (1. + val * val); // 1 + tan^2 = sec^2
        Dual::new(val, dval)
    }
}

// =============================================================================
// Exp & Logarithm
// =============================================================================

/// Exponential & Logarithm for Dual
impl ExpLogOps for Dual {
    type Output = Dual;

    fn exp(&self) -> Dual {
        let val = self.value().exp();
        let dval = val * self.slope();
        Dual::new(val, dval)
    }

    fn ln(&self) -> Dual {
        assert_ne!(self.value(), 0.);
        let val = self.value().ln();
        let dval = self.slope() / self.value();
        Dual::new(val, dval)
    }

    fn log(&self, base: f64) -> Self::Output {
        self.ln() / base.ln()
    }
}

// =============================================================================
// Power
// =============================================================================

/// Power for Dual
impl PowOps for Dual {
    type Output = Dual;

    fn pow(&self, n: usize) -> Dual {
        let x = self.value();
        let val = x.powi(n as i32);
        let dval = (n as f64) * x.powi((n - 1) as i32) * self.slope();
        Dual::new(val, dval)
    }

    fn powf(&self, f: f64) -> Dual {
        let x = self.value();
        let val = x.powf(f);
        let dval = f * x.powf(f - 1f64) * self.slope();
        Dual::new(val, dval)
    }

    fn sqrt(&self) -> Dual {
        self.powf(0.5)
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
        self.clone().into_iter().map(|x| f(x)).collect::<Vec<Dual>>()
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
        unimplemented!()
    }

    fn norm(&self) -> Self::Scalar {
        unimplemented!()
    }
}