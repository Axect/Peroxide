extern crate rand;
pub use self::rand::distributions::uniform::SampleUniform;
pub use self::rand::prelude::*;
pub use self::OPDist::*;
pub use self::TPDist::*;
use special::function::gaussian;
use statistics::rand::ziggurat;
use std::convert::Into;
use structure::matrix::*;
use structure::vector::*;


/// One parameter distribution
///
/// # Distributions
/// * `Bernoulli(prob)`: Bernoulli distribution
pub enum OPDist<T: PartialOrd + SampleUniform + Copy + Into<f64>> {
    Bernoulli(T)
}

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

/// RNG for OPDist
impl<T: PartialOrd + SampleUniform + Copy + Into<f64>> RNG for OPDist<T> {
    fn sample(&self, n: usize) -> Vec<f64> {
        match self {
            Bernoulli(prob) => {
                assert!((*prob).into() <= 1f64, "Probability should be smaller than 1");

                let mut rng = thread_rng();
                let mut v = vec![0f64; n];

                for i in 0..n {
                    let uniform = rng.gen_range(0f64, 1f64);
                    if uniform <= (*prob).into() {
                        v[i] = 1f64;
                    } else {
                        v[i] = 0f64;
                    }
                }
                v
            }
        }
    }

    fn pdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64 {
        match self {
            Bernoulli(prob) => {
                if x.into() == 1f64 {
                    (*prob).into()
                } else {
                    1f64 - (*prob).into()
                }
            }
        }
    }
}

/// RNG for TPDist
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
