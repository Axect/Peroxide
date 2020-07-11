use std::borrow;
pub use smallvec::Array;
use smallvec::{SmallVec, IntoIter};
use std::iter::FromIterator;
use std::slice::SliceIndex;
use std::ops::{Index, IndexMut, Deref, DerefMut, Neg, Add, Sub, Mul, Div};
use crate::statistics::ops::C;
use crate::traits::{
    fp::FPVector,
    num::{PowOps, ExpLogOps, TrigOps},
};

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct AD<A> where A: Array<Item=f64> {
    data: SmallVec<A>
}

// =============================================================================
// Inherit SmallVec
// =============================================================================
impl<A: Array<Item=f64> + Default> AD<A> {
    pub fn new() -> Self {
        AD {
            data: SmallVec::default()
        }
    }

    pub fn from_array(a: A) -> Self {
        AD {
            data: SmallVec::from(a)
        }
    }

    pub fn with_capacity(n: usize) -> Self {
        AD {
            data: SmallVec::with_capacity(n)
        }
    }

    pub fn set_len(&mut self, n: usize) {
        unsafe {
            self.data.set_len(n)
        }
    }

    pub fn empty(n: usize) -> Self {
        let mut x = Self::with_capacity(n);
        x.set_len(n);
        x
    }

    pub fn copy(x: &AD<A>) -> Self {
        let mut z = AD::empty(x.len());
        x.iter().enumerate().for_each(|(i, &t)| z[i] = t);
        z
    }
}

impl<A: Array<Item=f64> + Default> Default for AD<A> {
    fn default() -> Self {
        AD::new()
    }
}

unsafe impl<A: Array<Item=f64>> Array for AD<A> {
    type Item = f64;

    fn size() -> usize {
        A::size()
    }
}


impl<A: Array<Item=f64>> borrow::Borrow<[f64]> for AD<A> {
    fn borrow(&self) -> &[f64] {
        self.data.borrow()
    }
}

impl<A: Array<Item=f64>> borrow::BorrowMut<[f64]> for AD<A> {
    fn borrow_mut(&mut self) -> &mut [f64] {
        self.data.borrow_mut()
    }
}

impl<A: Array<Item=f64>> From<A> for AD<A> {
    fn from(x: A) -> Self {
        AD {
            data: SmallVec::from(x)
        }
    }
}

impl<A: Array<Item=f64> + Default> FromIterator<f64> for AD<A> {
    fn from_iter<T: IntoIterator<Item=f64>>(iter: T) -> Self {
        AD {
            data: SmallVec::from_iter(iter)
        }
    }
}

impl<A: Array<Item=f64>, I: SliceIndex<[f64]>> Index<I> for AD<A> {
    type Output = <I as SliceIndex<[f64]>>::Output;

    fn index(&self, index: I) -> &Self::Output {
        self.data.index(index)
    }
}

impl<A: Array<Item=f64>, I: SliceIndex<[f64]>> IndexMut<I> for AD<A> {
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        self.data.index_mut(index)
    }
}

impl<A: Array<Item=f64>> IntoIterator for AD<A> {
    type Item = f64;
    type IntoIter = IntoIter<A>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<A: Array<Item=f64>> Deref for AD<A> {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        self.data.deref()
    }
}

impl<A: Array<Item=f64>> DerefMut for AD<A> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.data.deref_mut()
    }
}

impl<A: Array<Item=f64> + Default + Clone> FPVector for AD<A> {
    type Scalar = f64;

    fn fmap<F>(&self, f: F) -> Self where
        F: Fn(Self::Scalar) -> Self::Scalar {
        self.iter().map(|&x| f(x)).collect()
    }

    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
        T: Into<Self::Scalar> {
        self.iter().fold(init.into(), |s, &x| f(s, x))
    }

    fn zip_with<F>(&self, f: F, other: &Self) -> Self where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar {
        self.iter().zip(other.iter()).map(|(&x, &y)| f(x, y)).collect()
    }

    fn filter<F>(&self, f: F) -> Self where
        F: Fn(Self::Scalar) -> bool {
        self.clone().into_iter().filter(|x| f(*x)).collect()
    }

    fn take(&self, n: usize) -> Self {
        self.clone().into_iter().take(n).collect()
    }

    fn skip(&self, n: usize) -> Self {
        self.clone().into_iter().skip(n).collect()
    }

    fn sum(&self) -> Self::Scalar {
        self.iter().sum()
    }

    fn prod(&self) -> Self::Scalar {
        self.iter().product()
    }
}

// =============================================================================
// Ops for AD
// =============================================================================
impl<A: Array<Item=f64> + Default> Neg for AD<A> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        self.iter().map(|x| -x).collect()
    }
}

impl<'a, 'b, A: Array<Item=f64> + Default> Add<&'b AD<A>> for &'a AD<A> {
    type Output = AD<A>;

    fn add(self, rhs: &'b AD<A>) -> Self::Output {
        self.data.iter().zip(rhs.data.iter()).map(|(x, y)| x + y).collect()
    }
}

impl<'a, 'b, A: Array<Item=f64> + Default> Sub<&'b AD<A>> for &'a AD<A> {
    type Output = AD<A>;

    fn sub(self, rhs: &'b AD<A>) -> Self::Output {
        self.iter().zip(rhs.iter()).map(|(x, y)| x - y).collect()
    }
}

impl<'a, 'b, A: Array<Item=f64> + Default> Mul<&'b AD<A>> for &'a AD<A> {
    type Output = AD<A>;

    fn mul(self, rhs: &'b AD<A>) -> Self::Output {
        mul(self, rhs)
    }
}

impl<'a, 'b, A: Array<Item=f64> + Default> Div<&'b AD<A>> for &'a AD<A> {
    type Output = AD<A>;

    fn div(self, rhs: &'b AD<A>) -> Self::Output {
        div(self, rhs)
    }
}

impl<A: Array<Item=f64> + Default> PowOps for AD<A> {
    fn powi(&self, n: i32) -> Self {
        powi(self, n as usize)
    }

    fn powf(&self, f: f64) -> Self {
        powf(self, f)
    }

    fn pow(&self, _f: Self) -> Self {
        unimplemented!()
    }

    fn sqrt(&self) -> Self {
        self.powf(0.5)
    }
}

impl<A: Array<Item=f64> + Default + Clone> ExpLogOps for AD<A> {
    fn exp(&self) -> Self {
        exp(self)
    }

    fn ln(&self) -> Self {
        ln(self)
    }

    fn log(&self, base: f64) -> Self {
        log(self, base)
    }

    fn log2(&self) -> Self {
        log(self, 2f64)
    }

    fn log10(&self) -> Self {
        log(self, 10f64)
    }
}

impl<A: Array<Item=f64> + Default> TrigOps for AD<A> {
    fn sin_cos(&self) -> (Self, Self) {
        sin_cos(self)
    }

    fn sin(&self) -> Self {
        let (s, _) = self.sin_cos();
        s
    }

    fn cos(&self) -> Self {
        let (_, c) = self.sin_cos();
        c
    }

    fn tan(&self) -> Self {
        tan(self)
    }

    fn sinh(&self) -> Self {
        let (s, _) = sinh_cosh(self);
        s
    }

    fn cosh(&self) -> Self {
        let (_, c) = sinh_cosh(self);
        c
    }

    fn tanh(&self) -> Self {
        tanh(self)
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

    fn asinh(&self) -> Self {
        unimplemented!()
    }

    fn acosh(&self) -> Self {
        unimplemented!()
    }

    fn atanh(&self) -> Self {
        unimplemented!()
    }
}

// =============================================================================
// Backend
// =============================================================================
fn mul<A: Array<Item=f64> + Default>(x: &AD<A>, y: &AD<A>) -> AD<A> {
    let mut z = AD::empty(x.len());
    for t in 0..z.len() {
        z[t] = if t < y.len() {
            x.iter()
                .take(t + 1)
                .zip(y.iter().take(t + 1).rev())
                .enumerate()
                .fold(0f64, |s, (i, (x1, y1))| s + (C(t, i) as f64) * x1 * y1)
        } else if t < x.len() {
            x.iter()
                .take(t + 1)
                .rev()
                .zip(y.iter())
                .enumerate()
                .fold(0f64, |s, (i, (x1, x2))| s + (C(t, i) as f64) * x1 * x2)
        } else {
            x.iter()
                .enumerate()
                .skip(t - y.len() + 1)
                .zip(y.iter().rev())
                .fold(0f64, |s, ((i, x1), y1)| s + (C(t, i) as f64) * x1 * y1)
        };
    }
    z
}

fn div<A: Array<Item=f64> + Default>(x: &AD<A>, y: &AD<A>) -> AD<A>{
    let mut z = AD::empty(x.len());
    z[0] = x[0] / y[0];
    let y0 = 1f64 / y[0];
    for i in 1..z.len() {
        let mut s = 0f64;
        for (j, (y1, z1)) in y[1..i + 1].iter().zip(z[0..i].iter().rev()).enumerate() {
            s += (C(i, j + 1) as f64) * y1 * z1;
        }
        z[i] = y0 * (x[i] - s);
    }
    z
}

fn ln<A: Array<Item=f64> + Default>(x: &AD<A>) -> AD<A> {
    let mut z = AD::empty(x.len());
    z[0] = x[0].ln();
    let x0 = 1f64 / x[0];
    for i in 1..z.len() {
        let mut s = 0f64;
        for (j, (x1, z1)) in x[1..i].iter().zip(z[1..i].iter().rev()).enumerate() {
            s += (C(i - 1, j + 1) as f64) * x1 * z1;
        }
        z[i] = x0 * (x[i] - s);
    }
    z
}

fn powf<A: Array<Item=f64> + Default>(x: &AD<A>, f: f64) -> AD<A> {
    let ln_x = ln(x);
    let mut z = AD::empty(x.len());
    z[0] = x[0].powf(f);
    for i in 1..z.len() {
        let mut s = 0f64;
        for (j, (z1, ln_x1)) in z[1..i].iter().zip(ln_x[1..i].iter().rev()).enumerate() {
            s += (C(i - 1, j + 1) as f64) * z1 * ln_x1;
        }
        z[i] = f * (z[0] * ln_x[i] + s);
    }
    z
}

fn powi<A: Array<Item=f64> + Default>(x: &AD<A>, n: usize) -> AD<A> {
    let mut z = AD::copy(x);
    for _i in 1..n {
        z = mul(&z, x);
    }
    z
}

fn exp<A: Array<Item=f64> + Default>(x: &AD<A>) -> AD<A> {
    let mut z = AD::empty(x.len());
    z[0] = x[0].exp();
    for i in 1..z.len() {
        z[i] = z[0..i]
            .iter()
            .zip(x[1..i + 1].iter().rev())
            .enumerate()
            .fold(0f64, |x, (k, (z1, x1))| x + (C(i - 1, k) as f64) * x1 * z1);
    }
    z
}

fn sin_cos<A: Array<Item=f64> + Default>(x: &AD<A>) -> (AD<A>, AD<A>) {
    let mut u = AD::empty(x.len());
    let mut v = AD::empty(x.len());
    u[0] = x[0].sin();
    v[0] = x[0].cos();

    for i in 1..x.len() {
        u[i] = v[0..i]
            .iter()
            .zip(x[1..i + 1].iter().rev())
            .enumerate()
            .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
        v[i] = -u[0..i]
            .iter()
            .zip(x[1..i + 1].iter().rev())
            .enumerate()
            .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
    }
    (u, v)
}

fn tan<A: Array<Item=f64> + Default>(x: &AD<A>) -> AD<A> {
    let (s, c) = sin_cos(x);
    div(&s, &c)
}

fn recip<A: Array<Item=f64> + Default>(x: &AD<A>) -> AD<A> {
    let mut z = AD::empty(x.len());
    z[0] = 1f64 / x[0];
    for i in 1..z.len() {
        let s = z[0..i]
            .iter()
            .zip(x[1..i + 1].iter().rev())
            .enumerate()
            .fold(0f64, |x, (k, (z1, x1))| x + (C(i, k) as f64) * z1 * x1);
        z[i] = -s * z[0];
    }
    z
}

#[allow(dead_code)]
fn csc_sec<A: Array<Item=f64> + Default>(x: &AD<A>) -> (AD<A>, AD<A>) {
    let (s, c) = sin_cos(x);
    (recip(&s), recip(&c))
}

#[allow(dead_code)]
fn cot<A: Array<Item=f64> + Default>(x: &AD<A>) -> AD<A> {
    let (s, c) = sin_cos(x);
    div(&c, &s)
}

fn log<A: Array<Item=f64> + Default + Clone>(x: &AD<A>, base: f64) -> AD<A> {
    ln(x).fmap(|t| t / base.ln())
}

pub fn powd<A: Array<Item=f64> + Default>(a: f64, x: &AD<A>) -> AD<A> {
    let mut z = AD::empty(x.len());
    z[0] = a.powf(x[0]);
    for i in 1..z.len() {
        z[i] = a.ln()
            * z[0..i]
                .iter()
                .zip(x[1..i + 1].iter().rev())
                .enumerate()
                .fold(0f64, |x, (k, (z1, x1))| x + (C(i - 1, k) as f64) * x1 * z1);
    }
    z
}

fn sinh_cosh<A: Array<Item=f64> + Default>(x: &AD<A>) -> (AD<A>, AD<A>) {
    let mut u = AD::empty(x.len());
    let mut v = AD::empty(x.len());
    u[0] = x[0].sinh();
    v[0] = x[0].cosh();

    for i in 1..x.len() {
        u[i] = v[0..i]
            .iter()
            .zip(x[1..i + 1].iter().rev())
            .enumerate()
            .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
        v[i] = u[0..i]
            .iter()
            .zip(x[1..i + 1].iter().rev())
            .enumerate()
            .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
    }
    (u, v)
}

fn tanh<A: Array<Item=f64> + Default>(x: &AD<A>) -> AD<A> {
    let (s, c) = sinh_cosh(x);
    div(&s, &c)
}

#[allow(dead_code)]
fn coth<A: Array<Item=f64> + Default>(x: &AD<A>) -> AD<A> {
    let (s, c) = sinh_cosh(x);
    div(&c, &s)
}

#[allow(dead_code)]
fn csch_sech<A: Array<Item=f64> + Default>(x: &AD<A>) -> (AD<A>, AD<A>) {
    let (s, c) = sinh_cosh(x);
    (recip(&s), recip(&c))
}
