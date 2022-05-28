//! Macro to non macro function

extern crate rand;
use self::rand::prelude::*;
use crate::structure::{
    matrix::Shape::{Col, Row},
    matrix::{matrix, Matrix, Shape},
};

use std::convert::Into;

/// R like seq function
///
/// # Example
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = seq(1, 10, 2);
/// assert_eq!(a, vec![1f64,3f64,5f64,7f64,9f64]);
/// 
/// let b = seq(1, 1, 1);
/// assert_eq!(b, vec![1f64]);
/// ```
pub fn seq<S, T, U>(start: S, end: T, step: U) -> Vec<f64>
where
    S: Into<f64> + Copy,
    T: Into<f64> + Copy,
    U: Into<f64> + Copy,
{
    let s = start.into();
    let e = end.into();
    let step = step.into();

    assert!(e >= s);

    let factor: f64 = (e - s) / step;
    let l: usize = factor as usize + 1;
    let mut v: Vec<f64> = vec![0f64; l];

    for i in 0..l {
        v[i] = s + step * (i as f64);
    }
    v
}

/// MATLAB like zeros (Matrix)
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = zeros(2, 2);
/// assert_eq!(a, matrix(vec![0f64;4], 2, 2, Row));
/// ```
pub fn zeros(r: usize, c: usize) -> Matrix {
    matrix(vec![0f64; r * c], r, c, Row)
}

/// Zeros with custom shape
pub fn zeros_shape(r: usize, c: usize, shape: Shape) -> Matrix {
    matrix(vec![0f64; r * c], r, c, shape)
}

pub fn concat<T: Clone + Copy>(v1: &Vec<T>, v2: &Vec<T>) -> Vec<T> {
    let mut v = v1.clone();
    v.extend_from_slice(&v2[..]);

    v
}

pub fn cat<T: Clone + Copy + Default>(val: T, vec: &Vec<T>) -> Vec<T> {
    let mut v = vec![val];
    v.extend_from_slice(&vec[..]);

    v
}

/// MATLAB like eye - Identity matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = eye(2);
/// assert_eq!(a, MATLAB::new("1 0;0 1"));
/// ```
pub fn eye(n: usize) -> Matrix {
    let mut m = zeros(n, n);
    for i in 0..n {
        m[(i, i)] = 1f64;
    }
    m
}

/// eye with custom shape
pub fn eye_shape(n: usize, shape: Shape) -> Matrix {
    let mut m = zeros_shape(n, n, shape);
    for i in 0..n {
        m[(i, i)] = 1f64;
    }
    m
}

/// MATLAB like linspace
/// 
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// 
/// fn main() {
///     let a = linspace(1, 10, 10);
///     assert_eq!(a, seq(1,10,1));
///     assert_eq!(a.len(), 10);
/// }
/// ```
pub fn linspace<S, T>(start: S, end: T, length: usize) -> Vec<f64>
where
    S: Into<f64> + Copy,
    T: Into<f64> + Copy,
{
    let step: f64 = (end.into() - start.into()) / (length as f64 - 1f64);
    seq(start, end, step)
}

/// MATLAB like logspace
/// 
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// 
/// fn main() {
///     let a = logspace(1, 1024, 11, 2);
///     let b = vec![1f64, 2f64, 4f64, 8f64, 16f64, 32f64, 64f64, 128f64, 256f64, 512f64, 1024f64];
///     assert_eq!(a, b);
/// }
/// ```
pub fn logspace<S, T, U>(start: S, end: T, length: usize, base: U) -> Vec<f64>
where
    S: Into<f64> + Copy,
    T: Into<f64> + Copy,
    U: Into<f64> + Copy,
{
    let s: f64 = start.into();
    let e: f64 = end.into();
    let b: f64 = base.into();

    assert!(e >= s);

    let step: f64 = (e.log(b) - s.log(b)) / (length as f64 - 1f64);

    let mut v: Vec<f64> = vec![0f64; length];

    for i in 0..length {
        v[i] = s * b.powf(step * (i as f64));
    }
    v
}

/// R like cbind - concatenate two matrix by column direction
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;4;1, 2, 2, Col);
///     let b = matrix!(5;8;1, 2, 2, Col);
///     let c = matrix!(1;8;1, 2, 4, Col);
///     assert_eq!(cbind(a,b), c);
/// }
/// ```
pub fn cbind(m1: Matrix, m2: Matrix) -> Matrix {
    let mut temp = m1;
    if temp.shape != Col {
        temp = temp.change_shape();
    }

    let mut temp2 = m2;
    if temp2.shape != Col {
        temp2 = temp2.change_shape();
    }

    let mut v = temp.data;
    let mut c = temp.col;
    let r = temp.row;

    assert_eq!(r, temp2.row);
    v.extend_from_slice(&temp2.data[..]);
    c += temp2.col;
    matrix(v, r, c, Col)
}

/// R like rbind - concatenate two matrix by row direction
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = matrix!(1;4;1, 2, 2, Row);
///     let b = matrix!(5;8;1, 2, 2, Row);
///     let c = matrix!(1;8;1, 4, 2, Row);
///     assert_eq!(rbind(a,b), c);
/// }
/// ```
pub fn rbind(m1: Matrix, m2: Matrix) -> Matrix {
    let mut temp = m1;
    if temp.shape != Row {
        temp = temp.change_shape();
    }

    let mut temp2 = m2;
    if temp2.shape != Row {
        temp2 = temp2.change_shape();
    }

    let mut v = temp.data;
    let c = temp.col;
    let mut r = temp.row;

    assert_eq!(c, temp2.col);
    v.extend_from_slice(&temp2.data[..]);
    r += temp2.row;
    matrix(v, r, c, Row)
}

/// Rand matrix
///
/// # Description
///
/// Range = from 0 to 1
pub fn rand(r: usize, c: usize) -> Matrix {
    let mut m = zeros(r, c);
    let mut rng = thread_rng();
    for i in 0..r {
        for j in 0..c {
            m[(i, j)] = rng.gen_range(0f64..=1f64);
        }
    }
    m
}
