//! Taylor mode forward automatic differentiation
//!
//! ## Important Features
//!
//! * Can automatic differentiate up to 2nd order (`AD0` ~ `AD2`)
//! * All `AD` are in stack (Guarantee high performance)
//! * You can see `AD` via `.print()`
//!
//! ## Implemented Traits for `AD`
//!
//! * `#[derive(Debug, Copy, Clone, PartialEq)]`
//! * `std::fmt::Display`
//! * `IntoIterator<Item = f64>`
//! * `IntoIterator<Item = &f64>`
//! * `IntoIterator<Item = &mut f64>`
//! * `FromIterator<f64>`
//! * `FromIterator<&f64>`
//! * `Index`, `IndexMut`
//! * `std::ops::{Neg, Add, Sub, Mul, Div}`
//! * `peroxide::traits::num::{PowOps, ExpLogOps, TrigOps}`
//! * `peroxide::util::print::Printable`
//!
//! ## Iterator of `AD`
//!
//! There are three iterators.
//!
//! * `ADIntoIter`
//! * `ADIter<'a>`
//! * `ADIterMut<'a>`
//!
//! Each implements `DoubleEndedIterator`, `ExactSizeIterator` also.
//!
//! ## Default Constructor
//!
//! * `AD0(f64)` : just constant
//! * `AD1(f64, f64)` : 1st order AD
//! * `AD2(f64, f64, f64)` : 2nd order AD
//!
//! ## Methods
//!
//! * `empty(&self) -> AD`
//! * `to_order(&self, n: usize) -> AD`
//! * `iter(&self) -> ADIter{i}`
//! * `iter_mut(&self) -> ADIterMut{i}`
//! * `len(&self) -> usize`
//! * `order(&self) -> usize`
//! * `from_order(n: usize) -> AD`
//! * `x(&self) -> f64`
//! * `dx(&self) -> f64`
//! * `ddx(&self) -> f64`
//!
//! ## Implemented Operations
//!
//! * `Add, Sub, Mul, Div`
//! * `sin, cos, tan`
//! * `sinh, cosh, tanh`
//! * `sin_cos`, `sinh_cosh`
//! * `exp, ln, log, log2, log10`
//! * `powi, powf, sqrt, pow`
//! * `asin`, `acos`, `atan`
//! * `asinh`, `acosh`, `atanh`
//!
//! ## Usage
//!
//! ### Construction
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // Declare x where x = 2
//!     let a = AD1(2f64, 1f64);
//!     // Declare x^2 where x = 2
//!     let b = AD2(4f64, 4f64, 2f64);
//!     // Convert AD1 -> AD2
//!     let c = a.to_order(2);
//!     // Zeros
//!     let d = AD::from_order(2);
//!
//!     assert_eq!(c, AD2(2f64, 1f64, 0f64));
//!     assert_eq!(d, AD2(0f64, 0f64, 0f64));
//! }
//! ```
//!
//! ### Operation
//!
//! For every binary operation, it returns higher order AD
//! (ex: AD1 + AD2 = AD2)
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = AD1(2f64, 1f64);       // x        at x = 2
//!     let b = AD2(4f64, 4f64, 2f64); // x^2      at x = 2
//!     let c = a + b;                 // x^2 + x  at x = 2
//!     let d = a * b;                 // x^3      at x = 2
//!     let e = a / b;                 // 1/x      at x = 2
//!     assert_eq!(c, AD2(6f64, 5f64, 2f64));
//!     assert_eq!(d, AD2(8f64, 12f64, 12f64));
//!     assert_eq!(e, AD2(0.5, -0.25, 0.25));
//! }
//! ```
//!

use crate::statistics::ops::C;
use crate::traits::{
    num::{ExpLogOps, PowOps, TrigOps},
    stable::StableFn,
    fp::FPVector,
    math::Vector,
    sugar::VecOps,
};
use std::iter::{FromIterator, DoubleEndedIterator, ExactSizeIterator};
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};
use self::AD::{AD0, AD1, AD2};

#[derive(Debug, Copy, Clone, PartialOrd, PartialEq)]
pub enum AD {
    AD0(f64),
    AD1(f64, f64),
    AD2(f64, f64, f64),
}

impl std::fmt::Display for AD {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!("{:?}", self);
        write!(f, "{}", s)
    }
}

impl AD {
    pub fn to_order(&self, n: usize) -> Self {
        if n == self.order() {
            return *self;
        }

        let mut z = match n {
            0 => AD0(0f64),
            1 => AD1(0f64, 0f64),
            2 => AD2(0f64, 0f64, 0f64),
            _ => panic!("No more index exists"),
        };

        for i in 0..z.len().min(self.len()) {
            z[i] = self[i];
        }

        z
    }

    pub fn order(&self) -> usize {
        match self {
            AD0(_) => 0,
            AD1(_, _) => 1,
            AD2(_, _, _) => 2,
        }
    }

    pub fn len(&self) -> usize {
        match self {
            AD0(_) => 1,
            AD1(_, _) => 2,
            AD2(_, _, _) => 3,
        }
    }

    pub fn iter(&self) -> ADIter {
        self.into_iter()
    }

    pub fn iter_mut(&mut self) -> ADIterMut {
        self.into_iter()
    }

    pub fn from_order(n: usize) -> Self {
        match n {
            0 => AD0(0f64),
            1 => AD1(0f64, 0f64),
            2 => AD2(0f64, 0f64, 0f64),
            _ => panic!("Not yet implemented higher order AD"),
        }
    }

    pub fn empty(&self) -> Self {
        match self {
            AD0(_) => AD0(0f64),
            AD1(_,_) => AD1(0f64, 0f64),
            AD2(_,_,_) => AD2(0f64, 0f64, 0f64),
        }
    }

    pub fn set_x(&mut self, x: f64) {
        match self {
            AD0(t) => {
                *t = x;
            }
            AD1(t, _) => {
                *t = x;
            }
            AD2(t, _, _) => {
                *t = x;
            }
        }
    }

    pub fn set_dx(&mut self, dx: f64) {
        match self {
            AD0(_) => panic!("Can't set dx for AD0"),
            AD1(_, dt) => {
                *dt = dx;
            }
            AD2(_, dt, _) => {
                *dt = dx;
            }
        }
    }

    pub fn set_ddx(&mut self, ddx: f64) {
        match self {
            AD0(_) => panic!("Can't set ddx for AD0"),
            AD1(_,_) => panic!("Can't set ddx for AD1"),
            AD2(_,_,ddt) => {
                *ddt = ddx;
            } 
        }
    }

    pub fn x(&self) -> f64 {
        match self {
            AD0(x) => *x,
            AD1(x, _) => *x,
            AD2(x, _, _) => *x,
        }
    }

    pub fn dx(&self) -> f64 {
        match self {
            AD0(_) => 0f64,
            AD1(_, dx) => *dx,
            AD2(_, dx, _) => *dx,
        }
    }

    pub fn ddx(&self) -> f64 {
        match self {
            AD0(_) => 0f64,
            AD1(_, _) => 0f64,
            AD2(_, _, ddx) => *ddx,
        }
    }

    pub fn x_ref(&self) -> Option<&f64> {
        match self {
            AD0(x) => Some(x),
            AD1(x, _) => Some(x),
            AD2(x, _, _) => Some(x),
        }
    }

    pub fn dx_ref(&self) -> Option<&f64> {
        match self {
            AD0(_) => None,
            AD1(_, dx) => Some(dx),
            AD2(_, dx, _) => Some(dx),
        }
    }

    pub fn ddx_ref(&self) -> Option<&f64> {
        match self {
            AD0(_) => None,
            AD1(_, _) => None,
            AD2(_, _, ddx) => Some(ddx),
        }
    }

    pub fn x_mut(&mut self) -> Option<&mut f64> {
        match self {
            AD0(x) => Some(x),
            AD1(x, _) => Some(x),
            AD2(x, _, _) => Some(x),
        }
    }

    pub fn dx_mut(&mut self) -> Option<&mut f64> {
        match self {
            AD0(_) => None,
            AD1(_, dx) => Some(dx),
            AD2(_, dx, _) => Some(dx),
        }
    }

    pub fn ddx_mut(&mut self) -> Option<&mut f64> {
        match self {
            AD0(_) => None,
            AD1(_, _) => None,
            AD2(_, _, ddx) => Some(ddx),
        }
    }

    #[allow(dead_code)]
    unsafe fn x_ptr(&self) -> Option<*const f64> {
        match self {
            AD0(x) => Some(&*x),
            AD1(x,_) => Some(&*x),
            AD2(x,_,_) => Some(&*x),
        }
    }

    #[allow(dead_code)]
    unsafe fn dx_ptr(&self) -> Option<*const f64> {
        match self {
            AD0(_) => None,
            AD1(_,dx) => Some(&*dx),
            AD2(_,dx,_) => Some(&*dx),
        }
    }

    #[allow(dead_code)]
    unsafe fn ddx_ptr(&self) -> Option<*const f64> {
        match self {
            AD0(_) => None,
            AD1(_,_) => None,
            AD2(_,_,ddx) => Some(&*ddx),
        }
    }

    unsafe fn x_mut_ptr(&mut self) -> Option<*mut f64> {
        match self {
            AD0(x) => Some(&mut *x),
            AD1(x,_) => Some(&mut *x),
            AD2(x,_,_) => Some(&mut *x),
        }
    }

    unsafe fn dx_mut_ptr(&mut self) -> Option<*mut f64> {
        match self {
            AD0(_) => None,
            AD1(_,dx) => Some(&mut *dx),
            AD2(_,dx,_) => Some(&mut *dx),
        }
    }

    unsafe fn ddx_mut_ptr(&mut self) -> Option<*mut f64> {
        match self {
            AD0(_) => None,
            AD1(_,_) => None,
            AD2(_,_,ddx) => Some(&mut *ddx),
        }
    }
}

impl Index<usize> for AD {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => self.x_ref().unwrap(),
            1 => self.dx_ref().expect("AD0 has no dx"),
            2 => self.ddx_ref().expect("AD0, AD1 have no ddx"),
            _ => panic!("No more index exists"),
        }
    }
}

impl IndexMut<usize> for AD {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => self.x_mut().unwrap(),
            1 => self.dx_mut().expect("AD0 has no dx"),
            2 => self.ddx_mut().expect("AD0, AD1 have no ddx"),
            _ => panic!("No more index exists"),
        }
    }
}

#[derive(Debug)]
pub struct ADIntoIter {
    ad: AD,
    index: usize,
    r_index: usize,
}

#[derive(Debug)]
pub struct ADIter<'a> {
    ad: &'a AD,
    index: usize,
    r_index: usize,
}

#[derive(Debug)]
pub struct ADIterMut<'a> {
    ad: &'a mut AD,
    index: usize,
    r_index: usize,
}

impl IntoIterator for AD {
    type Item = f64;
    type IntoIter = ADIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        ADIntoIter {
            ad: self,
            index: 0,
            r_index: 0,
        }
    }
}

impl<'a> IntoIterator for &'a AD {
    type Item = &'a f64;
    type IntoIter = ADIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        ADIter {
            ad: self,
            index: 0,
            r_index: 0,
        }
    }
}

impl<'a> IntoIterator for &'a mut AD {
    type Item = &'a mut f64;
    type IntoIter = ADIterMut<'a>;

    fn into_iter(self) -> Self::IntoIter {
        ADIterMut {
            ad: self,
            index: 0,
            r_index: 0,
        }
    }
}

impl Iterator for ADIntoIter {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.ad.len();
        if self.index + self.r_index < l {
            let result = match self.index {
                0 => Some(self.ad.x()),
                1 => match self.ad {
                    AD0(_) => None,
                    AD1(_, dx) => Some(dx),
                    AD2(_, dx, _) => Some(dx),
                },
                2 => match self.ad {
                    AD0(_) => None,
                    AD1(_, _) => None,
                    AD2(_, _, ddx) => Some(ddx),
                },
                _ => None,
            };
            self.index += 1;
            result
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let lower = self.ad.len() - (self.index + self.r_index);
        let upper = self.ad.len() - (self.index + self.r_index);
        (lower, Some(upper))
    }
}

impl<'a> Iterator for ADIter<'a> {
    type Item = &'a f64;

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.ad.len();
        if self.index + self.r_index < l {
            let result = match self.index {
                0 => self.ad.x_ref(),
                1 => self.ad.dx_ref(),
                2 => self.ad.ddx_ref(),
                _ => None,
            };
            self.index += 1;
            result
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let lower = self.ad.len() - (self.index + self.r_index);
        let upper = self.ad.len() - (self.index + self.r_index);
        (lower, Some(upper))
    }
}

impl<'a> Iterator for ADIterMut<'a> {
    type Item = &'a mut f64;

    fn next(&mut self) -> Option<Self::Item> {
        let l = self.ad.len();
        if self.index + self.r_index < l {
            unsafe {
                let result= match self.index {
                    0 => self.ad.x_mut_ptr(),
                    1 => self.ad.dx_mut_ptr(),
                    2 => self.ad.ddx_mut_ptr(),
                    _ => None,
                };
                self.index += 1;
                match result {
                    None => None,
                    Some(ad) => Some(&mut *ad)
                }
            }
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let lower = self.ad.len() - (self.index + self.r_index);
        let upper = self.ad.len() - (self.index + self.r_index);
        (lower, Some(upper))
    }
}

impl FromIterator<f64> for AD {
    fn from_iter<T: IntoIterator<Item = f64>>(iter: T) -> Self {
        let into_iter = iter.into_iter();
        let s = into_iter.size_hint().0 - 1;
        let mut z = match s {
            0 => AD0(0f64),
            1 => AD1(0f64, 0f64),
            2 => AD2(0f64, 0f64, 0f64),
            _ => panic!("Higher than order 3 is not allowed"),
        };
        for (i, elem) in into_iter.enumerate() {
            z[i] = elem;
        }
        z
    }
}

impl<'a> FromIterator<&'a f64> for AD {
    fn from_iter<T: IntoIterator<Item = &'a f64>>(iter: T) -> Self {
        let into_iter = iter.into_iter();
        let s = into_iter.size_hint().0 - 1;
        let mut z = match s {
            0 => AD0(0f64),
            1 => AD1(0f64, 0f64),
            2 => AD2(0f64, 0f64, 0f64),
            _ => panic!("Higher than order 3 is not allowed"),
        };
        for (i, &elem) in into_iter.enumerate() {
            z[i] = elem;
        }
        z
    }
}

impl DoubleEndedIterator for ADIntoIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.index + self.r_index == self.ad.len() {
            return None;
        }
        let order = self.ad.order();
        let result = self.ad[order - self.r_index];
        self.r_index += 1;
        Some(result)
    }
}

impl<'a> DoubleEndedIterator for ADIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.index + self.r_index == self.ad.len() {
            return None;
        }
        let order = self.ad.order();
        let result = &self.ad[order - self.r_index];
        self.r_index += 1;
        Some(result)
    }
}

impl ExactSizeIterator for ADIntoIter {
    fn len(&self) -> usize {
        self.ad.len() - (self.index + self.r_index)
    }
}

impl<'a> ExactSizeIterator for ADIter<'a> {
    fn len(&self) -> usize {
        self.ad.len() - (self.index + self.r_index)
    }
}

impl Neg for AD {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.into_iter().map(|x| -x).collect()
    }
}

impl Add<AD> for AD {
    type Output = Self;

    fn add(self, rhs: AD) -> Self::Output {
        let ord = self.order().max(rhs.order());
        let (a, b) = (self.to_order(ord), rhs.to_order(ord));

        a.into_iter()
            .zip(b.into_iter())
            .map(|(x, y)| x + y)
            .collect()
    }
}

impl Sub<AD> for AD {
    type Output = Self;

    fn sub(self, rhs: AD) -> Self::Output {
        let ord = self.order().max(rhs.order());
        let (a, b) = (self.to_order(ord), rhs.to_order(ord));

        a.into_iter()
            .zip(b.into_iter())
            .map(|(x, y)| x - y)
            .collect()
    }
}

impl Mul<AD> for AD {
    type Output = Self;

    fn mul(self, rhs: AD) -> Self::Output {
        let ord = self.order().max(rhs.order());
        let (a, b) = (self.to_order(ord), rhs.to_order(ord));

        let mut z = a.clone();
        for t in 0..z.len() {
            z[t] = a
                .into_iter()
                .take(t + 1)
                .zip(b.into_iter().take(t + 1).rev())
                .enumerate()
                .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
        }
        z
    }
}

impl Div<AD> for AD {
    type Output = Self;

    fn div(self, rhs: AD) -> Self::Output {
        let ord = self.order().max(rhs.order());
        let (a, b) = (self.to_order(ord), rhs.to_order(ord));

        let mut z = a.clone();
        z[0] = a[0] / b[0];
        let y0 = 1f64 / b[0];
        for i in 1 .. z.len() {
            let mut s = 0f64;
            for (j, (&y1, &z1)) in b.iter().skip(1).take(i).zip(z.iter().take(i).rev()).enumerate() {
                s += (C(i, j+1) as f64) * y1 * z1;
            }
            z[i] = y0 * (a[i] - s);
        }
        z
    }
}

impl ExpLogOps for AD {
    fn exp(&self) -> Self {
        let mut z = self.empty();
        z[0] = self[0].exp();
        for i in 1 .. z.len() {
            z[i] = z.iter()
                .take(i)
                .zip(self.iter().skip(1).take(i).rev())
                .enumerate()
                .fold(0f64, |x, (k, (&z1, &x1))| x + (C(i-1, k) as f64) * x1 * z1);
        }
        z
    }

    fn ln(&self) -> Self {
        let mut z = self.empty();
        z[0] = self[0].ln();
        let x0 = 1f64 / self[0];
        for i in 1 .. z.len() {
            let mut s = 0f64;
            for (k, (&z1, &x1)) in z.iter().skip(1).take(i-1).zip(self.iter().skip(1).take(i-1).rev()).enumerate() {
                s += (C(i-1, k+1) as f64) * z1 * x1;
            }
            z[i] = x0 * (self[i] - s);
        }
        z
    }

    fn log(&self, base: f64) -> Self {
        self.ln().iter().map(|x| x / base.ln()).collect()
    }
}

impl PowOps for AD {
    fn powi(&self, n: i32) -> Self {
        let mut z = self.clone();
        for _i in 1 .. n {
            z = z * *self;
        }
        z
    }

    fn powf(&self, f: f64) -> Self {
        let ln_x = self.ln();
        let mut z = self.empty();
        z[0] = self.x().powf(f);
        for i in 1 .. z.len() {
            let mut s = 0f64;
            for (j, (&z1, &ln_x1)) in z.iter().skip(1).take(i-1).zip(ln_x.iter().skip(1).take(i-1).rev()).enumerate() {
                s += (C(i-1, j+1) as f64) * z1 * ln_x1;
            }
            z[i] = f * (z[0] * ln_x[i] + s);
        }
        z
    }

    fn pow(&self, y: Self) -> Self {
        let ln_x = self.ln();
        let p = y * ln_x;
        let mut z = self.empty();
        z[0] = self.x().powf(y.x());
        for n in 1 .. z.len() {
            let mut s = 0f64;
            for (k, (&z1, &p1)) in z.iter().skip(1).take(n-1).zip(p.iter().skip(1).take(n-1).rev()).enumerate() {
                s += (C(n-1, k+1) as f64) * z1 * p1; 
            }
            z[n] = z[0] * p[n] + s;
        }
        z
    }
}

impl TrigOps for AD {
    fn sin_cos(&self) -> (Self, Self) {
        let mut u = self.empty();
        let mut v = self.empty();
        u[0] = self[0].sin();
        v[0] = self[0].cos();
        for i in 1 .. u.len() {
            u[i] = v.iter()
                .take(i)
                .zip(self.iter().skip(1).take(i).rev())
                .enumerate()
                .fold(0f64, |x, (k, (&v1, &x1))| x + (C(i-1, k) as f64) * x1 * v1);
            v[i] = u.iter()
                .take(i)
                .zip(self.iter().skip(1).take(i).rev())
                .enumerate()
                .fold(0f64, |x, (k, (&u1, &x1))| x + (C(i-1, k) as f64) * x1 * u1);
        }
        (u, v)
    }

    fn sinh_cosh(&self) -> (Self, Self) {
        let mut u = self.empty();
        let mut v = self.empty();
        u[0] = self[0].sinh();
        v[0] = self[0].cosh();
        for i in 1 .. u.len() {
            u[i] = v.iter()
                .take(i)
                .zip(self.iter().skip(1).take(i).rev())
                .enumerate()
                .fold(0f64, |x, (k, (&v1, &x1))| x + (C(i-1, k) as f64) * x1 * v1);
            v[i] = u.iter()
                .take(i)
                .zip(self.iter().skip(1).take(i).rev())
                .enumerate()
                .fold(0f64, |x, (k, (&u1, &x1))| x + (C(i-1, k) as f64) * x1 * u1);
        }
        (u, v)
    }

    fn asin(&self) -> Self {
        let dx = 1f64 / (1f64 - self.powi(2)).sqrt();
        let mut z = self.empty();
        z[0] = self[0].asin();
        for n in 1 .. z.len() {
            z[n] = dx.iter()
                .take(n)
                .zip(self.iter().skip(1).take(n).rev())
                .enumerate()
                .fold(0f64, |s, (k, (&q1, &x1))| s + (C(n-1, k) as f64) * x1 * q1);
        }
        z
    }

    fn acos(&self) -> Self {
        let dx = (-1f64) / (1f64 - self.powi(2)).sqrt();
        let mut z = self.empty();
        z[0] = self[0].acos();
        for n in 1 .. z.len() {
            z[n] = dx.iter()
                .take(n)
                .zip(self.iter().skip(1).take(n).rev())
                .enumerate()
                .fold(0f64, |s, (k, (&q1, &x1))| s + (C(n-1, k) as f64) * x1 * q1);
        }
        z
    }

    fn atan(&self) -> Self {
        let dx = 1f64 / (1f64 + self.powi(2));
        let mut z = self.empty();
        z[0] = self[0].atan();
        for n in 1 .. z.len() {
            z[n] = dx.iter()
                .take(n)
                .zip(self.iter().skip(1).take(n).rev())
                .enumerate()
                .fold(0f64, |s, (k, (&q1, &x1))| s + (C(n-1, k) as f64) * x1 * q1);
        }
        z
    }

    fn asinh(&self) -> Self {
        let dx = 1f64 / (1f64 + self.powi(2)).sqrt();
        let mut z = self.empty();
        z[0] = self[0].asinh();
        for n in 1 .. z.len() {
            z[n] = dx.iter()
                .take(n)
                .zip(self.iter().skip(1).take(n).rev())
                .enumerate()
                .fold(0f64, |s, (k, (&q1, &x1))| s + (C(n-1, k) as f64) * x1 * q1);
        }
        z
    }

    fn acosh(&self) -> Self {
        let dx = 1f64 / (self.powi(2) - 1f64).sqrt();
        let mut z = self.empty();
        z[0] = self[0].acosh();
        for n in 1 .. z.len() {
            z[n] = dx.iter()
                .take(n)
                .zip(self.iter().skip(1).take(n).rev())
                .enumerate()
                .fold(0f64, |s, (k, (&q1, &x1))| s + (C(n-1, k) as f64) * x1 * q1);
        }
        z
    }

    fn atanh(&self) -> Self {
        let dx = 1f64 / (1f64 - self.powi(2));
        let mut z = self.empty();
        z[0] = self[0].atanh();
        for n in 1 .. z.len() {
            z[n] = dx.iter()
                .take(n)
                .zip(self.iter().skip(1).take(n).rev())
                .enumerate()
                .fold(0f64, |s, (k, (&q1, &x1))| s + (C(n-1, k) as f64) * x1 * q1);
        }
        z
    }
}

impl From<f64> for AD {
    fn from(other: f64) -> Self {
        AD0(other)
    }
}

impl From<AD> for f64 {
    fn from(other: AD) -> Self {
        other.x()
    }
}

impl Add<f64> for AD {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        let mut z = self;
        z[0] += rhs;
        z
    }
}

impl Sub<f64> for AD {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        let mut z = self;
        z[0] -= rhs;
        z
    }
}

impl Mul<f64> for AD {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        self.iter().map(|&x| x * rhs).collect()
    }
}

impl Div<f64> for AD {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        self.iter().map(|&x| x / rhs).collect()
    }
}

impl Add<AD> for f64 {
    type Output = AD;

    fn add(self, rhs: AD) -> Self::Output {
        let mut z = rhs;
        z[0] += self;
        z
    }
}

impl Sub<AD> for f64 {
    type Output = AD;

    fn sub(self, rhs: AD) -> Self::Output {
        let mut z = rhs.empty();
        z[0] = self;
        z - rhs
    }
}

impl Mul<AD> for f64 {
    type Output = AD;

    fn mul(self, rhs: AD) -> Self::Output {
        rhs.iter().map(|&x| x * self).collect()
    }
}

impl Div<AD> for f64 {
    type Output = AD;

    fn div(self, rhs: AD) -> Self::Output {
        let ad0 = AD::from(self);
        ad0 / rhs
    }
}

// ad_struct_def!();
// ad_display!();
// ad_impl!();
// ad_impl_from!();
// ad_iter_def!();
// ad_impl_into_iter!();
// ad_impl_iter!();
// ad_impl_from_iter!();
// ad_impl_double_ended_iter!();
// ad_impl_exact_size_iter!();
// ad_impl_index!();
// ad_impl_neg!();
// ad_impl_add!();
// ad_impl_sub!();
// ad_impl_mul!();
// ad_impl_div!();
// ad_impl_explogops!();
// ad_impl_powops!();
// ad_impl_trigops!();
// ad_impl_ad!();
// def_ad!();
// ad_impl_from_type!(f64);
// ad_impl_from_type!(f32);
// ad_impl_from_type!(i64);
// ad_impl_from_type!(i32);
// ad_impl_add_f64!();
// ad_impl_sub_f64!();
// ad_impl_mul_f64!();
// ad_impl_div_f64!();
// f64_impl_add_ad!();
// f64_impl_sub_ad!();
// f64_impl_mul_ad!();
// f64_impl_div_ad!();
// f64_impl_from_ad!();
// ad_impl_stable_fn!();

// pub trait AD:
//     std::fmt::Display
//     + Clone
//     + Copy
//     + PartialEq
//     + From<AD1>
//     + From<AD2>
//     + From<AD3>
//     + From<AD4>
//     + From<AD5>
//     + From<AD6>
//     + From<AD7>
//     + From<AD8>
//     + From<AD9>
//     + From<AD10>
//     + Into<AD1>
//     + Into<AD2>
//     + Into<AD3>
//     + Into<AD4>
//     + Into<AD5>
//     + Into<AD6>
//     + Into<AD7>
//     + Into<AD8>
//     + Into<AD9>
//     + Into<AD10>
//     + From<f64>
//     + Into<f64>
//     + Add<Output = Self>
//     + Sub<Output = Self>
//     + Mul<Output = Self>
//     + Div<Output = Self>
//     + Add<f64, Output = Self>
//     + Sub<f64, Output = Self>
//     + Mul<f64, Output = Self>
//     + Div<f64, Output = Self>
//     + PowOps
//     + ExpLogOps
//     + TrigOps
// {
//     fn to_ad1(self) -> AD1 {
//         self.into()
//     }
//
//     fn to_ad2(self) -> AD2 {
//         self.into()
//     }
//
//     fn to_ad3(self) -> AD3 {
//         self.into()
//     }
//
//     fn to_ad4(self) -> AD4 {
//         self.into()
//     }
//
//     fn to_ad5(self) -> AD5 {
//         self.into()
//     }
//
//     fn to_ad6(self) -> AD6 {
//         self.into()
//     }
//
//     fn to_ad7(self) -> AD7 {
//         self.into()
//     }
//
//     fn to_ad8(self) -> AD8 {
//         self.into()
//     }
//
//     fn to_ad9(self) -> AD9 {
//         self.into()
//     }
//
//     fn to_ad10(self) -> AD10 {
//         self.into()
//     }
// }

//impl AD for f64 {}

/// Generic AD functions
///
/// # Description
/// To lift `AD` functions
///
/// # Implementation
///
/// * All `Fn(AD) -> AD` functions can be lift to `Fn(f64) -> f64` via `StableFn<f64>`
/// * `grad(&self) -> Self` gives gradient of original function
/// * But still can use `Fn(AD) -> AD` via `StableFn<AD>`
///
/// # Usage
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let ad0 = 2f64;
///     let ad1 = AD1(2f64, 1f64);
///     let ad2 = AD2(2f64, 1f64, 0f64);
///
///     let f_ad = ADFn::new(f);
///
///     let ans0 = ad0.powi(2);
///     let ans1 = ad1.powi(2).dx();
///     let ans2 = ad2.powi(2).ddx();
///
///     assert_eq!(ans0, f_ad.call_stable(ad0));
///
///     let df = f_ad.grad();
///     assert_eq!(ans1, df.call_stable(ad0));
///
///     let ddf = df.grad();
///     assert_eq!(ans2, ddf.call_stable(ad0));
///
///     let ad1_output = f_ad.call_stable(ad1);
///     assert_eq!(ad1_output, AD1(4f64, 4f64));
///
///     let ad2_output = f_ad.call_stable(ad2);
///     assert_eq!(ad2_output, AD2(4f64, 4f64, 2f64));
/// }
///
/// fn f(x: AD) -> AD {
///     x.powi(2)
/// }
/// ```
pub struct ADFn<F> {
    f: Box<F>,
    grad_level: usize,
}

impl<F: Clone> ADFn<F> {
    pub fn new(f: F) -> Self {
        Self {
            f: Box::new(f),
            grad_level: 0usize,
        }
    }

    /// Gradient
    pub fn grad(&self) -> Self {
        assert!(self.grad_level < 2, "Higher order AD is not allowed");
        ADFn {
            f: (self.f).clone(),
            grad_level: self.grad_level + 1,
        }
    }
}

impl<F: Fn(AD) -> AD> StableFn<f64> for ADFn<F> {
    type Output = f64;
    fn call_stable(&self, target: f64) -> Self::Output {
        match self.grad_level {
            0 => (self.f)(AD::from(target)).into(),
            1 => (self.f)(AD1(target, 1f64)).dx(),
            2 => (self.f)(AD2(target, 1f64, 0f64)).ddx(),
            _ => panic!("Higher order AD is not allowed"),
        }
    }
}

impl<F: Fn(AD) -> AD> StableFn<AD> for ADFn<F> {
    type Output = AD;
    fn call_stable(&self, target: AD) -> Self::Output {
        (self.f)(target)
    }
}

impl<F: Fn(Vec<AD>) -> Vec<AD>> StableFn<Vec<f64>> for ADFn<F> {
    type Output = Vec<f64>;
    fn call_stable(&self, target: Vec<f64>) -> Self::Output {
        ((self.f)(target.iter().map(|&t| AD::from(t)).collect())).iter().map(|&t| t.x()).collect()
    }
}

impl<F: Fn(Vec<AD>) -> Vec<AD>> StableFn<Vec<AD>> for ADFn<F> {
    type Output = Vec<AD>;
    fn call_stable(&self, target: Vec<AD>) -> Self::Output {
        (self.f)(target)
    }
}

// Nightly only
//pub struct ADLift<F>(F);
//
//impl<F: FnOnce<(T)>, T> FnOnce<f64> for ADLift<F> {
//    type Output = f64;
//
//}

pub trait ADVec {
    fn to_ad_vec(&self) -> Vec<AD>;
    fn to_f64_vec(&self) -> Vec<f64>;
}

impl ADVec for Vec<f64> {
    fn to_ad_vec(&self) -> Vec<AD> {
        self.iter().map(|&t| AD::from(t)).collect()
    }

    fn to_f64_vec(&self) -> Vec<f64> {
        self.iter().map(|&t| t).collect()
    }
}

impl ADVec for Vec<AD> {
    fn to_ad_vec(&self) -> Vec<AD> {
        self.iter().map(|&t| t).collect()
    }

    fn to_f64_vec(&self) -> Vec<f64> {
        self.iter().map(|t| t.x()).collect()
    }
}

impl FPVector for Vec<AD> {
    type Scalar = AD;

    fn fmap<F>(&self, f: F) -> Self
    where
            F: Fn(Self::Scalar) -> Self::Scalar {
        self.iter().map(|&x| f(x)).collect()
    }

    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
            F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
            T: Into<Self::Scalar> {
        self.iter().fold(init.into(), |x, &y| f(x,y))
    }

    fn filter<F>(&self, f: F) -> Self
    where
            F: Fn(Self::Scalar) -> bool {
        self.iter().filter(|&x| f(*x)).map(|&t| t).collect()
    }

    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
            F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar {
        self.iter().zip(other.iter()).map(|(&x, &y)| f(x, y)).collect()
    }

    fn take(&self, n: usize) -> Self {
        self.iter().take(n).map(|&t| t).collect()
    }

    fn skip(&self, n: usize) -> Self {
        self.iter().skip(n).map(|&t| t).collect()
    }

    fn sum(&self) -> Self::Scalar {
        let s = self[0];
        self.reduce(s, |x, y| x + y)
    }

    fn prod(&self) -> Self::Scalar {
        let s = self[0];
        self.reduce(s, |x, y| x * y)
    }
}

impl Vector for Vec<AD> {
    fn add_vec<'a, 'b>(&'a self, rhs: &'b Self) -> Self {
        self.add_v(rhs)
    }

    fn sub_vec<'a, 'b>(&'a self, rhs: &'b Self) -> Self {
        self.sub_v(rhs)
    }

    fn mul_scalar<T: Into<f64>>(&self, rhs: T) -> Self {
        self.mul_s(AD::from(rhs.into()))
    }
}

impl VecOps for Vec<AD> {}
