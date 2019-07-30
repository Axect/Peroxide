use operation::number::Number::{D, E, F};
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::process::exit;
use {Dual, Real};
use {ExpLogOps, TrigOps};
use {HyperDual, PowOps};

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum NumError {
    DiffType,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Number {
    F(f64),
    D(Dual),
    E(NumError),
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
            E(x) => E(x),
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
            (E(x), _) => E(x),
            (_, E(y)) => E(y),
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
            (E(x), _) => E(x),
            (_, E(y)) => E(y),
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
            (E(x), _) => E(x),
            (_, E(y)) => E(y),
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
            (E(x), _) => E(x),
            (_, E(y)) => E(y),
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
            E(x) => E(x.to_owned()),
        }
    }

    fn ln(&self) -> Self {
        match self {
            F(x) => F(x.ln()),
            D(x) => D(x.exp()),
            E(x) => E(x.to_owned()),
        }
    }

    fn log(&self, base: f64) -> Self {
        match self {
            F(x) => F(x.log(base)),
            D(x) => D(x.log(base)),
            E(x) => E(x.to_owned()),
        }
    }

    fn log2(&self) -> Self {
        match self {
            F(x) => F(x.log2()),
            D(x) => D(x.log2()),
            E(x) => E(x.to_owned()),
        }
    }

    fn log10(&self) -> Self {
        match self {
            F(x) => F(x.log10()),
            D(x) => D(x.log10()),
            E(x) => E(x.to_owned()),
        }
    }
}

impl TrigOps for Number {
    fn sin(&self) -> Self {
        match self {
            F(x) => F(x.sin()),
            D(x) => D(x.sin()),
            E(x) => E(x.to_owned()),
        }
    }

    fn cos(&self) -> Self {
        match self {
            F(x) => F(x.cos()),
            D(x) => D(x.cos()),
            E(x) => E(x.to_owned()),
        }
    }

    fn tan(&self) -> Self {
        match self {
            F(x) => F(x.tan()),
            D(x) => D(x.tan()),
            E(x) => E(x.to_owned()),
        }
    }

    fn asin(&self) -> Self {
        match self {
            F(x) => F(x.asin()),
            D(x) => D(x.asin()),
            E(x) => E(x.to_owned()),
        }
    }

    fn acos(&self) -> Self {
        match self {
            F(x) => F(x.acos()),
            D(x) => D(x.acos()),
            E(x) => E(x.to_owned()),
        }
    }

    fn atan(&self) -> Self {
        match self {
            F(x) => F(x.atan()),
            D(x) => D(x.atan()),
            E(x) => E(x.to_owned()),
        }
    }

    fn sinh(&self) -> Self {
        match self {
            F(x) => F(x.sinh()),
            D(x) => D(x.sinh()),
            E(x) => E(x.to_owned()),
        }
    }

    fn cosh(&self) -> Self {
        match self {
            F(x) => F(x.cosh()),
            D(x) => D(x.cosh()),
            E(x) => E(x.to_owned()),
        }
    }

    fn tanh(&self) -> Self {
        match self {
            F(x) => F(x.tanh()),
            D(x) => D(x.tanh()),
            E(x) => E(x.to_owned()),
        }
    }

    fn asinh(&self) -> Self {
        match self {
            F(x) => F(x.asinh()),
            D(x) => D(x.asinh()),
            E(x) => E(x.to_owned()),
        }
    }

    fn acosh(&self) -> Self {
        match self {
            F(x) => F(x.acosh()),
            D(x) => D(x.acosh()),
            E(x) => E(x.to_owned()),
        }
    }

    fn atanh(&self) -> Self {
        match self {
            F(x) => F(x.atanh()),
            D(x) => D(x.atanh()),
            E(x) => E(x.to_owned()),
        }
    }

    fn sin_cos(&self) -> (Self, Self) {
        match self {
            F(x) => (F(x.sin()), F(x.cos())),
            D(x) => (D(x.sin()), D(x.cos())),
            E(x) => (E(x.to_owned()), E(x.to_owned())),
        }
    }
}

impl PowOps for Number {
    fn powi(&self, n: i32) -> Self {
        match self {
            F(x) => F(x.powi(n)),
            D(x) => D(x.powi(n)),
            E(x) => E(x.to_owned()),
        }
    }

    #[allow(unreachable_patterns)]
    fn powf(&self, f: Self) -> Self {
        match (self, &f) {
            (F(x), F(y)) => F(x.powf(*y)),
            (F(x), D(y)) => D(Dual::from_f64(*x).powf(*y)),
            (D(x), F(y)) => D(x.powf(Dual::from_f64(*y))),
            (D(x), D(y)) => D(x.powf(*y)),
            (E(x), _) => E(x.to_owned()),
            (_, E(y)) => E(y.to_owned()),
            _ => unreachable!()
        }
    }

    fn sqrt(&self) -> Self {
        match self {
            F(x) => F(x.sqrt()),
            D(x) => D(x.sqrt()),
            E(x) => E(x.to_owned()),
        }
    }
}

impl Real for Number {
    fn to_f64(&self) -> f64 {
        match self {
            F(x) => x.to_owned(),
            D(x) => x.to_f64(),
            E(x) => {
                eprintln!("error {:?}", x.to_owned());
                exit(1);
            }
        }
    }

    fn from_f64(f: f64) -> Self {
        F(f)
    }

    fn to_dual(&self) -> Dual {
        match self {
            F(x) => x.to_dual(),
            D(x) => x.to_owned(),
            E(x) => {
                eprintln!("error {:?}", x.to_owned());
                exit(1);
            }
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
            E(x) => E(x),
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
        self.clone().into_iter().map(|x| x.to_hyper_dual()).collect()
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