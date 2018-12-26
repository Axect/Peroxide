#[allow(unused_imports)]
use structure::matrix::*;
#[allow(unused_imports)]
use structure::vector::*;
#[allow(unused_imports)]
use statistics::stat::*;

use std::convert::Into;

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

pub fn zeros(r: usize, c: usize) -> Matrix {
    matrix(vec![0f64; r * c], r, c, Row)
}