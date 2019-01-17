#[allow(unused_imports)]
use structure::matrix::*;
#[allow(unused_imports)]
use structure::vector::*;
#[allow(unused_imports)]
use statistics::stat::*;

use std::convert::Into;

/// R like seq function
///
/// # Example
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = seq(1, 10, 2);
/// assert_eq!(a, vec![1f64,3f64,5f64,7f64,9f64]);
/// ```
pub fn seq<S,T,U>(start: S, end: T, step: U) -> Vector
    where S: Into<f64> + Copy,
          T: Into<f64> + Copy,
          U: Into<f64> + Copy {
    let s = start.into();
    let e = end.into();
    let step = step.into();

    assert!(e > s);

    let factor: f64 = (e - s) / step;
    let l: usize = factor as usize + 1;
    let mut v: Vec<f64> = Vec::new();
    
    for i in 0 .. l {
        v.push(s + step * (i as f64));
    }
    v
}

/// MATLAB like zeros (Matrix)
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = zeros(2, 2);
/// assert_eq!(a, matrix(vec![0f64;4], 2, 2, Row));
/// ```
pub fn zeros(r: usize, c: usize) -> Matrix {
    matrix(vec![0f64; r * c], r, c, Row)
}

pub fn concat<T: Clone + Copy + Default>(v1: Vec<T>, v2: Vec<T>) -> Vec<T> {
    let l1 = v1.len();
    let l2 = v2.len();

    let mut v = vec![Default::default(); l1 + l2];

    for i in 0 .. l1 {
        v[i] = v1[i];
    }

    for i in 0 .. l2 {
        v[i + l1] = v2[i];
    }
    v
}

/// MATLAB like eye - Identity matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = eye(2);
/// assert_eq!(a, MATLAB::new("1 0;0 1"));
/// ```
pub fn eye(n: usize) -> Matrix {
    let mut m = zeros(n, n);
    for i in 0 .. n {
        m[(i, i)] = 1f64;
    }
    m
}

/// R like cbind - concatenate two matrix by column direction
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Col);
/// let b = matrix!(5;8;1, 2, 2, Col);
/// let c = matrix!(1;8;1, 2, 4, Col);
/// assert_eq!(cbind(a,b), c);
/// ```
pub fn cbind(m1: Matrix, m2: Matrix) -> Matrix {
    let mut temp = m1;
    if temp.shape != Col {
        temp = temp.change_shape();
    }

    let mut v = temp.data;
    let mut c = temp.col;
    let r = temp.row;

    assert_eq!(r, m2.row);
    v.extend(&m2.data.clone());
    c += m2.col;
    matrix(v, r, c, Col)
}