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
    Normal,
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

pub trait RNG {
    fn sample(&self, n: usize) -> Vec<f64>;
}

impl<T> RNG for Rand<T> where T: convert::Into<f64> + SampleUniform + PartialOrd + Copy {
    fn sample(&self, n: usize) -> Vec<f64> {
        let mut rng = thread_rng();
        let mut v = vec![0f64; n];
        
        let start = self.range.0;
        let end = self.range.1;

        match self.dist {
            Uniform => {
                for i in 0 .. n {
                    v[i] = rng.gen_range(start, end).into();
                }
            },
            Normal => {
                for i in 0 .. n {
                    v[i] = rng.gen_range(start, end).into();
                }
            }
        }
        v
    } 
}