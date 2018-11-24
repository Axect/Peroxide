/// Structure for Automatic Differentiation

use std::fmt;
use std::convert::Into;
use std::ops::{Neg, Add, Sub, Mul, Div};
use operation::extra_ops::{TrigOps, ExpLogOps, PowOps};

/// Dual - Structure for AD
///
/// # Fields
/// `x: f64` : value
/// `dx: f64` : slope at `x`
///
/// But they are private fields.
/// You should use `value` or `slope` function to extract them.
#[derive(Debug, Copy, Clone)]
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
        Dual { x: x.into(), dx: dx.into() }
    }

    /// Just extract value
    pub fn value(&self) -> f64 {
        self.x
    }

    /// Just extract slope
    pub fn slope(&self) -> f64 {
        self.dx
    }
}

pub fn dual<T: Into<f64> + Copy>(x: T, dx: T) -> Dual {
    Dual::new(x, dx)
}

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
        assert_eq!(other.x, 0f64);
        let v1 = self.x;
        let v2 = other.x;
        let dv1 = self.dx;
        let dv2 = other.dx;

        Dual::new(v1 / v2, (dv1 * v2 - v1 * dv2) / (v2 * v2))
    }
}

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
}

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
        let dval = f * x.powf(f) * self.slope();
        Dual::new(val, dval)
    }
}