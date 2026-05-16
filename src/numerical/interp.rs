#[allow(unused_imports)]
use crate::structure::polynomial::*;

use std::convert::Into;
use std::f64::consts::PI;

pub fn chebyshev_nodes<T>(num: usize, start: T, end: T) -> Vec<f64>
where
    T: Into<f64> + Copy,
{
    let mut v = vec![0f64; num];
    let a = start.into();
    let b = end.into();
    for (i, slot) in v.iter_mut().enumerate() {
        *slot = (a + b) / 2. + 0.5 * (b - a) * ((2 * i + 1) as f64 * PI / (2 * num) as f64).cos();
    }
    v
}
