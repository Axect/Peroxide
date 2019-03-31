use operation::extra_ops::{ExpLogOps, PowOps, TrigOps};
use std::convert::Into;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};
use structure::vector::*;
use structure::dual::dual;

#[derive(Debug, Copy, Clone, Default)]
pub struct HyperDual {
    x: f64,
    dx: f64,
    ddx: f64,
}

impl fmt::Display for HyperDual {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s1 = format!(
            "value: {}\nslope: {}\nacc: {}",
            self.x,
            self.dx,
            self.ddx
        );
        write!(f, "{}", s1)
    }
}

impl HyperDual {
    pub fn new<T: Into<f64> + Copy>(x: T, dx: T, ddx: T) -> HyperDual {
        HyperDual {
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
    type Output = HyperDual;
    fn neg(self) -> HyperDual {
        HyperDual::new(-self.x, -self.dx, -self.ddx)
    }
}

impl Add<HyperDual> for HyperDual {
    type Output = HyperDual;

    fn add(self, rhs: HyperDual) -> Self::Output {
        HyperDual::new(
            self.x + rhs.x,
            self.dx + rhs.dx,
            self.ddx + rhs.ddx,
        )
    }
}

impl Sub<HyperDual> for HyperDual {
    type Output = HyperDual;

    fn sub(self, rhs: HyperDual) -> Self::Output {
        HyperDual::new(
            self.x - rhs.x,
            self.dx - rhs.dx,
            self.ddx - rhs.ddx,
        )
    }
}

impl Mul<HyperDual> for HyperDual {
    type Output = HyperDual;

    fn mul(self, rhs: HyperDual) -> Self::Output {
        let (x, dx, ddx) = self.extract();
        let (y, dy, ddy) = rhs.extract();

        HyperDual::new(
            x * y,
            dx*y + x*dy,
            ddx*y + 2f64*dx*dy + x*ddy
        )
    }
}

impl Div<HyperDual> for HyperDual {
    type Output = HyperDual;

    fn div(self, rhs: HyperDual) -> Self::Output {
        assert_ne!(rhs.x, 0f64);
        let (x, dx, ddx) = self.extract();
        let (y, dy, ddy) = rhs.extract();

        let dual_x = dual(x, dx);
        let dual_y = dual(y, dy);

        let x_div_y = (dual_x / dual_y).slope();

        HyperDual::new(
            x/y,
            (dx*y - x*dy)/(y*y),
            (ddx - 2f64*x_div_y*dy - x/y*ddy)/y
        )
    }
}

impl Add<f64> for HyperDual {
    type Output = HyperDual;

    fn add(self, rhs: f64) -> Self::Output {
        self + HyperDual::new(rhs, 0., 0.)
    }
}

impl Sub<f64> for HyperDual {
    type Output = HyperDual;

    fn sub(self, rhs: f64) -> Self::Output {
        self - HyperDual::new(rhs, 0., 0.)
    }
}

impl Mul<f64> for HyperDual {
    type Output = HyperDual;

    fn mul(self, rhs: f64) -> Self::Output {
        self * HyperDual::new(rhs, 0., 0.)
    }
}

impl Div<f64> for HyperDual {
    type Output = HyperDual;

    fn div(self, rhs: f64) -> Self::Output {
        self / HyperDual::new(rhs, 0., 0.)
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
    type Output = HyperDual;

    fn sin(&self) -> Self::Output {
        let x = self.x.sin();
        let dx = self.dx * self.x.cos();
        let ddx = self.ddx * self.x.cos() - self.dx.powi(2) * self.x.sin();
        HyperDual::new(x, dx, ddx)
    }

    fn cos(&self) -> Self::Output {
        let x = self.x.cos();
        let dx = - self.dx * self.x.sin();
        let ddx =  - self.ddx * self.x.sin() - self.dx.powi(2) * self.x.cos();
        HyperDual::new(x, dx, ddx)
    }

    fn tan(&self) -> Self::Output {
        unimplemented!()
    }
}