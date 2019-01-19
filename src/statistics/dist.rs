extern crate rand;
pub use self::rand::distributions::uniform::SampleUniform;
pub use self::rand::prelude::*;
pub use self::TPDist::*;
use special::function::gaussian;
use statistics::rand::ziggurat;
use std::convert::Into;
use structure::matrix::*;
use structure::vector::*;

/// Two parameter distribution
///
/// # Distributions
/// * `Uniform(start, end)`: Uniform distribution
/// * `Normal(mean, std)`: Normal distribution
pub enum TPDist<T: PartialOrd + SampleUniform + Copy + Into<f64>> {
    Uniform(T, T),
    Normal(T, T),
}

/// Random Number Generator trait
///
/// # Methods
/// * `sample`: extract samples
pub trait RNG {
    /// Extract samples of distributions
    fn sample(&self, n: usize) -> Vec<f64>;

    /// Probability Distribution Function
    ///
    /// # Type
    /// `f64 -> f64`
    fn pdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64;
}

impl<T: PartialOrd + SampleUniform + Copy + Into<f64>> RNG for TPDist<T> {
    fn sample(&self, n: usize) -> Vec<f64> {
        match self {
            Uniform(start, end) => {
                let mut rng = thread_rng();
                let mut v = vec![0f64; n];

                for i in 0..n {
                    v[i] = rng.gen_range(*start, *end).into();
                }
                v
            }
            Normal(m, s) => {
                let mut rng = thread_rng();
                let mut v = vec![0f64; n];

                for i in 0..n {
                    v[i] = ziggurat(&mut rng, (*s).into()) + (*m).into();
                }
                v
            }
        }
    }

    fn pdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64 {
        match self {
            Uniform(a, b) => {
                let length = (*b).into() - (*a).into();
                1f64 / length
            }
            Normal(m, s) => {
                let mean = (*m).into();
                let std = (*s).into();
                gaussian(x.into(), mean, std)
            }
        }
    }
}
