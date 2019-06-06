use std::ops::{Add, Mul, Div, Sub};

pub trait PowOps: Sized {
    fn powi(&self, n: i32) -> Self;
    fn powf(&self, f: f64) -> Self;
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

pub trait Real: PowOps + TrigOps + ExpLogOps + Add<Output=Self> + Mul<Output=Self> + Div<Output=Self> + Sub<Output=Self> {
    fn to_f64(&self) -> f64;
    fn from_f64(f: f64) -> Self;
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
}