//! Hyper dual number system for Automatic differentiation
//!
//! ## Hyper dual number system
//!
//! * `HyperDual` is structure for 2nd order AD
//!     * `value(&self) -> f64`
//!     * `slope(&self) -> f64`
//!     * `accel(&self) -> f64`
//! * There are two constructors
//!     * `HyperDual::new(T, T, T)`
//!     * `hyper_dual(T, T, T)`
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let x = hyper_dual(1, 1, 0); // x at x = 1
//!         (x.clone() + x.clone()).print(); // hyper_dual(2, 2, 0)
//!         (x.clone() * x.clone()).print(); // hyper_dual(1, 2, 2)
//!         // and etc.
//!     }
//!     ```
//!
//!     * Also, after `0.10.1`, you can use reference ops.
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         let x = hyper_dual(1,1,0);
//!         (&x + &x).print();
//!         (&x * &x).print();
//!         // and etc.
//!     }
//!     ```

use operation::extra_ops::{ExpLogOps, PowOps, TrigOps};
use std::convert::Into;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};
use structure::dual::dual;
#[allow(unused_imports)]
use structure::vector::*;

/// Hyper Dual number
///
/// # Description
///
/// For second order differentiation
#[derive(Debug, Copy, Clone, Default)]
pub struct HyperDual {
    x: f64,
    dx: f64,
    ddx: f64,
}

impl fmt::Display for HyperDual {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s1 = format!("value: {}\nslope: {}\naccel: {}", self.x, self.dx, self.ddx);
        write!(f, "{}", s1)
    }
}

impl HyperDual {
    pub fn new<T: Into<f64> + Copy>(x: T, dx: T, ddx: T) -> Self {
        Self {
            x: x.into(),
            dx: dx.into(),
            ddx: ddx.into(),
        }
    }

    pub fn value(&self) -> f64 {
        self.x
    }

    pub fn slope(&self) -> f64 {
        self.dx
    }

    pub fn accel(&self) -> f64 {
        self.ddx
    }

    pub fn extract(&self) -> (f64, f64, f64) {
        (self.x, self.dx, self.ddx)
    }
}

pub fn hyper_dual<T: Into<f64> + Copy>(x: T, dx: T, ddx: T) -> HyperDual {
    HyperDual::new(x, dx, ddx)
}

impl Neg for HyperDual {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.dx, -self.ddx)
    }
}

impl Add<HyperDual> for HyperDual {
    type Output = Self;

    fn add(self, rhs: HyperDual) -> Self::Output {
        Self::new(self.x + rhs.x, self.dx + rhs.dx, self.ddx + rhs.ddx)
    }
}

impl Sub<HyperDual> for HyperDual {
    type Output = Self;

    fn sub(self, rhs: HyperDual) -> Self::Output {
        Self::new(self.x - rhs.x, self.dx - rhs.dx, self.ddx - rhs.ddx)
    }
}

impl Mul<HyperDual> for HyperDual {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let (x, dx, ddx) = self.extract();
        let (y, dy, ddy) = rhs.extract();

        Self::new(x * y, dx * y + x * dy, ddx * y + 2f64 * dx * dy + x * ddy)
    }
}

impl Div<HyperDual> for HyperDual {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        assert_ne!(rhs.x, 0f64);
        let (x, dx, ddx) = self.extract();
        let (y, dy, ddy) = rhs.extract();

        let dual_x = dual(x, dx);
        let dual_y = dual(y, dy);

        let x_div_y = (dual_x / dual_y).slope();

        Self::new(
            x / y,
            (dx * y - x * dy) / (y * y),
            (ddx - 2f64 * x_div_y * dy - x / y * ddy) / y,
        )
    }
}

impl<'a, 'b> Add<&'b HyperDual> for &'a HyperDual {
    type Output = HyperDual;

    fn add(self, rhs: &HyperDual) -> Self::Output {
        HyperDual {
            x: self.x + rhs.x,
            dx: self.dx + rhs.dx,
            ddx: self.ddx + rhs.ddx,
        }
    }
}

impl<'a, 'b> Sub<&'b HyperDual> for &'a HyperDual {
    type Output = HyperDual;

    fn sub(self, rhs: &HyperDual) -> Self::Output {
        HyperDual {
            x: self.x - rhs.x,
            dx: self.dx - rhs.dx,
            ddx: self.ddx - rhs.ddx,
        }
    }
}

impl<'a, 'b> Mul<&'b HyperDual> for &'a HyperDual {
    type Output = HyperDual;

    fn mul(self, rhs: &HyperDual) -> Self::Output {
        HyperDual {
            x: self.x * rhs.x,
            dx: self.dx * rhs.x + self.x * rhs.dx,
            ddx: self.ddx * rhs.x + 2f64 * self.dx * rhs.dx + self.x * rhs.ddx,
        }
    }
}

impl<'a, 'b> Div<&'b HyperDual> for &'a HyperDual {
    type Output = HyperDual;

    fn div(self, rhs: &HyperDual) -> Self::Output {
        let (x, dx, ddx) = self.extract();
        let (y, dy, ddy) = rhs.extract();

        let dual_x = dual(x, dx);
        let dual_y = dual(y, dy);

        let x_div_y = (dual_x / dual_y).slope();

        HyperDual::new(
            x / y,
            (dx * y - x * dy) / (y * y),
            (ddx - 2f64 * x_div_y * dy - x / y * ddy) / y,
        )
    }
}

impl Add<f64> for HyperDual {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        self + Self::new(rhs, 0., 0.)
    }
}

impl Sub<f64> for HyperDual {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        self - Self::new(rhs, 0., 0.)
    }
}

impl Mul<f64> for HyperDual {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        self * Self::new(rhs, 0., 0.)
    }
}

impl Div<f64> for HyperDual {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        self / Self::new(rhs, 0., 0.)
    }
}

impl<'a> Add<f64> for &'a HyperDual {
    type Output = HyperDual;

    fn add(self, rhs: f64) -> Self::Output {
        self.add(&HyperDual::new(rhs, 0f64, 0f64))
    }
}

impl<'a> Sub<f64> for &'a HyperDual {
    type Output = HyperDual;

    fn sub(self, rhs: f64) -> Self::Output {
        self.sub(&HyperDual::new(rhs, 0f64, 0f64))
    }
}

impl<'a> Mul<f64> for &'a HyperDual {
    type Output = HyperDual;

    fn mul(self, rhs: f64) -> Self::Output {
        self.mul(&HyperDual::new(rhs, 0f64, 0f64))
    }
}

impl<'a> Div<f64> for &'a HyperDual {
    type Output = HyperDual;

    fn div(self, rhs: f64) -> Self::Output {
        self.div(&HyperDual::new(rhs, 0f64, 0f64))
    }
}

impl Add<HyperDual> for f64 {
    type Output = HyperDual;

    fn add(self, rhs: HyperDual) -> Self::Output {
        rhs.add(self)
    }
}

impl Sub<HyperDual> for f64 {
    type Output = HyperDual;

    fn sub(self, rhs: HyperDual) -> Self::Output {
        -rhs.sub(self)
    }
}

impl Mul<HyperDual> for f64 {
    type Output = HyperDual;

    fn mul(self, rhs: HyperDual) -> Self::Output {
        rhs.mul(self)
    }
}

impl Div<HyperDual> for f64 {
    type Output = HyperDual;

    fn div(self, rhs: HyperDual) -> Self::Output {
        hyper_dual(self, 0., 0.) / rhs
    }
}

impl TrigOps for HyperDual {
    fn sin(&self) -> Self {
        let x = self.x.sin();
        let dx = self.dx * self.x.cos();
        let ddx = self.ddx * self.x.cos() - self.dx.powi(2) * x;
        Self::new(x, dx, ddx)
    }

    fn cos(&self) -> Self {
        let x = self.x.cos();
        let dx = -self.dx * self.x.sin();
        let ddx = -self.ddx * self.x.sin() - self.dx.powi(2) * x;
        Self::new(x, dx, ddx)
    }

    fn tan(&self) -> Self {
        let x = self.x.tan();
        let dx = self.dx * (1f64 + x.powi(2));
        let ddx = self.ddx * (1f64 + x.powi(2)) + dx * self.dx * 2f64 * x;
        Self::new(x, dx, ddx)
    }

    fn asin(&self) -> Self {
        unimplemented!()
    }

    fn acos(&self) -> Self {
        unimplemented!()
    }

    fn atan(&self) -> Self {
        unimplemented!()
    }

    fn sinh(&self) -> Self {
        unimplemented!()
    }

    fn cosh(&self) -> Self {
        unimplemented!()
    }

    fn tanh(&self) -> Self {
        unimplemented!()
    }

    fn asinh(&self) -> Self {
        unimplemented!()
    }

    fn acosh(&self) -> Self {
        unimplemented!()
    }

    fn atanh(&self) -> Self {
        unimplemented!()
    }

    fn sin_cos(&self) -> (Self, Self) {
        unimplemented!()
    }
}

impl ExpLogOps for HyperDual {
    fn exp(&self) -> Self {
        let x = self.x.exp();
        let dx = self.dx * x;
        let ddx = self.ddx * x + self.dx.powi(2) * x;
        Self::new(x, dx, ddx)
    }

    fn ln(&self) -> Self {
        assert!(self.x > 0f64, "Logarithm Domain Error");
        let x = self.x.ln();
        let dx = self.dx / self.x;
        let ddx = self.ddx / self.x - self.dx.powi(2) / self.x.powi(2);
        Self::new(x, dx, ddx)
    }

    fn log(&self, base: f64) -> Self {
        self.ln() / base.ln()
    }

    fn log2(&self) -> Self {
        unimplemented!()
    }

    fn log10(&self) -> Self {
        unimplemented!()
    }
}

impl PowOps for HyperDual {
    fn powi(&self, n: i32) -> Self {
        let mut s = self.clone();
        for _i in 1..n {
            s = s * s;
        }
        s
    }

    fn powf(&self, f: Self) -> Self {
        let n = f.value();
        let dn = f.slope();
        let ddn = f.accel();

        if dn == 0f64 && ddn == 0f64 {
            let x = self.x.powf(n);
            let dx = self.dx * n * self.x.powf(n - 1f64);
            let ddx = self.ddx * n * self.x.powf(n - 1f64)
                + self.dx.powi(2) * n * (n - 1f64) * self.x.powf(n - 2f64);
            Self::new(x, dx, ddx)
        } else {
            unimplemented!()
        }
    }

    fn sqrt(&self) -> Self {
        self.powf(HyperDual::new(0.5f64, 0f64, 0f64))
    }
}
