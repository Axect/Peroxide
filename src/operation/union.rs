use ::{Dual, Real};
use std::ops::{Add, Sub, Mul};
use operation::union::Number::{D, F, E};
use operation::union::NumError::DiffType;

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum NumError {
    DiffType
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Number {
    F(f64),
    D(Dual),
    E(NumError)
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

impl<'a, 'b> Add<&'b Number> for &'a Number {
    type Output = Number;

    fn add(self, rhs: &Number) -> Self::Output {
        match (self, rhs) {
            (F(x), F(y)) => F(*x + *y),
            (D(x), D(y)) => D(x + y),
            (F(x), D(y)) => D(*x + *y),
            (D(x), F(y)) => D(*x + *y),
            (E(x), _) => E(x.to_owned()),
            (_, E(y)) => E(y.to_owned()),
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

impl<'a> Add<f64> for &'a Number {
    type Output = Number;

    fn add(self, rhs: f64) -> Self::Output {
        self.add(&F(rhs))
    }
}

impl<'a> Add<Dual> for &'a Number {
    type Output = Number;

    fn add(self, rhs: Dual) -> Self::Output {
        self.add(&D(rhs))
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

impl<'a, 'b> Sub<&'b Number> for &'a Number {
    type Output = Number;

    fn sub(self, rhs: &Number) -> Self::Output {
        match (self, rhs) {
            (F(x), F(y)) => F(*x - *y),
            (D(x), D(y)) => D(x - y),
            (F(x), D(y)) => D(*x - *y),
            (D(x), F(y)) => D(*x - *y),
            (E(x), _) => E(x.to_owned()),
            (_, E(y)) => E(y.to_owned())
        }
    }
}

impl Sub<f64> for Number {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        self.sub(F(rhs))
    }
}

impl<'a> Sub<f64> for &'a Number {
    type Output = Number;

    fn sub(self, rhs: f64) -> Self::Output {
        self.sub(&F(rhs))
    }
}

impl Sub<Dual> for Number {
    type Output = Self;

    fn sub(self, rhs: Dual) -> Self::Output {
        self.sub(D(rhs))
    }
}

impl<'a> Sub<Dual> for &'a Number {
    type Output = Number;

    fn sub(self, rhs: Dual) -> Self::Output {
        self.sub(&D(rhs))
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

impl<'a, 'b> Mul<&'b Number> for &'a Number {
    type Output = Number;

    fn mul(self, rhs: &Number) -> Self::Output {
        match (self, rhs) {
            (F(x), F(y)) => F(*x * *y),
            (D(x), D(y)) => D(x * y),
            (F(x), D(y)) => D(*x * *y),
            (D(x), F(y)) => D(*x * *y),
            (E(x), _) => E(x.to_owned()),
            (_, E(y)) => E(y.to_owned()),
        }
    }
}

impl Mul<f64> for Number {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        self.mul(F(rhs))
    }
}

impl<'a> Mul<f64> for &'a Number {
    type Output = Number;

    fn mul(self, rhs: f64) -> Self::Output {
        self.mul(&F(rhs))
    }
}

impl Mul<Dual> for Number {
    type Output = Self;

    fn mul(self, rhs: Dual) -> Self::Output {
        self.mul(D(rhs))
    }
}

impl<'a> Mul<Dual> for &'a Number {
    type Output = Number;

    fn mul(self, rhs: Dual) -> Self::Output {
        self.mul(&D(rhs))
    }
}