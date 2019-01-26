extern crate rand;
pub use self::rand::distributions::uniform::SampleUniform;
pub use self::rand::prelude::*;
pub use self::OPDist::*;
pub use self::TPDist::*;
use special::function::*;
use statistics::rand::ziggurat;
use statistics::stat::Statistics;
use std::convert::Into;
use structure::matrix::*;
use structure::vector::*;

/// One parameter distribution
///
/// # Distributions
/// * `Bernoulli(prob)`: Bernoulli distribution
pub enum OPDist<T: PartialOrd + SampleUniform + Copy + Into<f64>> {
    Bernoulli(T),
}

/// Two parameter distribution
///
/// # Distributions
/// * `Uniform(start, end)`: Uniform distribution
/// * `Normal(mean, std)`: Normal distribution
pub enum TPDist<T: PartialOrd + SampleUniform + Copy + Into<f64>> {
    Uniform(T, T),
    Normal(T, T),
    Beta(T, T),
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
                assert!(
                    (*prob).into() <= 1f64,
                    "Probability should be smaller than 1"
                );

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
            Beta(a, b) => {
                let mut rng1 = thread_rng();
                let mut rng2 = thread_rng();
                let mut v = vec![0f64; n];

                let a_f64 = (*a).into();
                let b_f64 = (*b).into();

                // For acceptance-rejection method
                let c_x = (a_f64 - 1f64) / (a_f64 + b_f64 - 2f64);
                let c = self.pdf(c_x); // Beta(mode(x) | a, b)

                let mut iter_num = 0usize;

                while iter_num < n {
                    let u1 = rng1.gen_range(0f64, 1f64);
                    let u2 = rng2.gen_range(0f64, 1f64);

                    if u2
                        <= 1f64 / (c * beta(a_f64, b_f64))
                            * u1.powf(a_f64 - 1f64)
                            * (1f64 - u1).powf(b_f64 - 1f64)
                    {
                        v[iter_num] = u1;
                        iter_num += 1;
                    }
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
            Beta(a, b) => {
                let a_f64 = (*a).into();
                let b_f64 = (*b).into();
                1f64 / beta(a_f64, b_f64)
                    * x.into().powf(a_f64 - 1f64)
                    * (1f64 - x.into()).powf(b_f64 - 1f64)
            }
        }
    }
}

impl<T: PartialOrd + SampleUniform + Copy + Into<f64>> Statistics for TPDist<T> {
    type Array = Vec<f64>;
    type Value = f64;

    fn mean(&self) -> Self::Value {
        match self {
            Uniform(a, b) => ((*a).into() + (*b).into()) / 2f64,
            Normal(m, _s) => (*m).into(),
            Beta(a, b) => (*a).into() / ((*a).into() + (*b).into()),
        }
    }

    fn var(&self) -> Self::Value {
        match self {
            Uniform(a, b) => ((*b).into() - (*a).into()).powi(2) / 12f64,
            Normal(_m, s) => (*s).into().powi(2),
            Beta(a, b) => {
                let a_f64 = (*a).into();
                let b_f64 = (*b).into();
                a_f64 * b_f64 / ((a_f64 + b_f64).powi(2) * (a_f64 + b_f64 + 1f64))
            }
        }
    }

    fn sd(&self) -> Self::Value {
        match self {
            Uniform(a, b) => self.var().sqrt(),
            Normal(_m, s) => (*s).into(),
            Beta(a, b) => self.var().sqrt(),
        }
    }

    fn cov(&self) -> Self::Array {
        unimplemented!()
    }

    fn cor(&self) -> Self::Array {
        unimplemented!()
    }
}