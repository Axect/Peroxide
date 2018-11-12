extern crate rand;
pub use self::rand::prelude::*;
pub use self::rand::distributions::uniform::SampleUniform;

use std::convert;
pub use self::Dist::*;

#[allow(unused_imports)]
use matrix::*;

#[derive(Debug, Clone, Copy)]
pub enum Dist {
    Uniform,
}

#[derive(Debug, Clone)]
pub struct Rand<T> where T: PartialOrd + SampleUniform + Copy{
    pub range: (T, T),
    pub dist: Dist,
}

impl<T> Rand<T> where T: PartialOrd + SampleUniform + Copy {
    pub fn new(range: (T, T), dist: Dist) -> Rand<T> {
        Rand {
            range: range,
            dist: dist,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Normal {
    pub m: f64,
    pub s: f64,
}

impl Normal {
    pub fn new_norm(m: f64, s: f64) -> Normal {
        Normal { m: m, s: s }
    }
}

pub trait RNG {
    fn sample(&self, n: usize) -> Vec<f64>;
}

impl<T> RNG for Rand<T> where T: convert::Into<f64> + SampleUniform + PartialOrd + Copy {
    fn sample(&self, n: usize) -> Vec<f64> {
        let mut rng = thread_rng();
        let mut v = vec![0f64; n];

        match self.dist {
            Uniform => {
                for i in 0 .. n {
                    let start = self.range.0;
                    let end = self.range.1;
                    v[i] = rng.gen_range(start, end).into();
                }
            },
        }
        v
    } 
}

impl RNG for Normal {
    fn sample(&self, n: usize) -> Vec<f64> {
        let mut rng = thread_rng();
        let mut v = vec![0f64; n];

        for i in 0 .. n {
            v[i] = marsaglia_polar(&mut rng, self.m, self.s);
        }
        v
    }
}

// =============================================================================
// Back end utils
// =============================================================================
pub fn marsaglia_polar(rng: &mut ThreadRng, m: f64, s: f64) -> f64 {
    let mut x1 = 0f64;
    let mut x2 = 0f64;
    let mut _y2 = 0f64;
    let mut w = 0f64;

    while w == 0. || w >= 1. {
        x1 = 2.0 * rng.gen_range(0f64, 1f64) - 1.0;
        x2 = 2.0 * rng.gen_range(0f64, 1f64) - 1.0;
        w = x1 * x1 + x2 * x2;
    }

    w = (-2.0 * w.ln() / w).sqrt();
    let y1 = x1 * w;
    _y2 = x2 * w;

    return m + y1 * s;
}