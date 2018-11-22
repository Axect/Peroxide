use std::fmt;
use std::convert::Into;
use std::ops::{Neg, Add, Sub, Mul, Div};
use operation::extra_ops::TrigOps;

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
    pub fn new<T>(x: T, dx: T) -> Dual where T: Into<f64> + Copy {
        Dual { x: x.into(), dx: dx.into() }
    }

    pub fn value(&self) -> f64 {
        self.x
    }

    pub fn slope(&self) -> f64 {
        self.dx
    }
}

impl Neg for Dual {
    type Output = Dual;
    fn neg(self) -> Dual {
        Dual::new(-self.x, -self.dx)
    }
}

impl Add<Dual> for Dual {
    type Output = Dual;
    fn add(self, other: Dual) -> Dual {
        Dual::new(self.x + other.x, self.dx + other.dx)
    }
}

impl Sub<Dual> for Dual {
    type Output = Dual;
    fn sub(self, other: Dual) -> Dual {
        Dual::new(self.x - other.x, self.dx - other.dx)
    }
}

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