//! Probabilistic distributions
//!
//! ## Probability Distribution
//!
//! * There are some famous pdf in Peroxide (not checked pdfs will be implemented soon)
//!     * Bernoulli
//!     * Beta
//!     * Dirichlet
//!     * Gamma
//!     * Normal
//!     * Student's t
//!     * Uniform
//! * There are two enums to represent probability distribution
//!     * `OPDist<T>` : One parameter distribution (Bernoulli)
//!     * `TPDist<T>` : Two parameter distribution (Uniform, Normal, Beta, Gamma)
//!         * `T: PartialOrd + SampleUniform + Copy + Into<f64>`
//! * There are some traits for pdf
//!     * `RNG` trait - extract sample & calculate pdf
//!     * `Statistics` trait - already shown above
//!
//! ### `RNG` trait
//!
//! * `RNG` trait is composed of two fields
//!     * `sample`: Extract samples
//!     * `pdf` : Calculate pdf value at specific point
//!     ```rust
//!     extern crate rand;
//!     use rand::distributions::uniform::SampleUniform;
//!
//!     pub trait RNG {
//!         /// Extract samples of distributions
//!         fn sample(&self, n: usize) -> Vec<f64>;
//!
//!         /// Probability Distribution Function
//!         ///
//!         /// # Type
//!         /// `f64 -> f64`
//!         fn pdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64;
//!     }
//!     ```
//!
//! ### Bernoulli Distribution
//!
//! * Definition
//!     $$ \text{Bern}(x | \mu) = \mu^x (1-\mu)^{1-x} $$
//! * Representative value
//!     * Mean: $\mu$
//!     * Var : $\mu(1 - \mu)$
//! * In peroxide, to generate $\text{Bern}(x | \mu)$, use simple traits
//!     1. Generate $U \sim \text{Unif}(0, 1)$
//!     2. If $U \leq \mu$, then $X = 1$ else $X = 0$
//! * Usage is very simple
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let b = Bernoulli(0.1); // Bern(x | 0.1)
//!         b.sample(100).print();  // Generate 100 samples
//!         b.pdf(0).print();       // 0.9
//!         b.mean().print();       // 0.1
//!         b.var().print();        // 0.09 (approximately)
//!         b.sd().print();         // 0.3  (approximately)
//!     }
//!     ```
//!
//! ### Uniform Distribution
//!
//! * Definition
//!     $$\text{Unif}(x | a, b) = \begin{cases} \frac{1}{b - a} & x \in [a,b]\\\ 0 & \text{otherwise} \end{cases}$$
//! * Representative value
//!     * Mean: $\frac{a + b}{2}$
//!     * Var : $\frac{1}{12}(b-a)^2$
//! * To generate uniform random number, Peroxide uses `rand` crate
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         // Uniform(start, end)
//!         let a = Uniform(0, 1);
//!         a.sample(100).print();
//!         a.pdf(0.2).print();
//!         a.mean().print();
//!         a.var().print();
//!         a.sd().print();
//!     }
//!     ```
//!
//! ### Normal Distribution
//!
//! * Definition
//!     $$\mathcal{N}(x | \mu, \sigma^2) = \frac{1}{\sqrt{2\pi \sigma^2}} \exp{\left( - \frac{(x - \mu)^2}{2\sigma^2}\right)}$$
//! * Representative value
//!     * Mean: $\mu$
//!     * Var: $\sigma^2$
//! * To generate normal random number, there are two famous algorithms
//!     * Marsaglia-Polar method
//!     * Ziggurat traits
//! * In peroxide (after ver 0.19.1), use `rand_distr` to generate random normal samples.
//! * <del>In peroxide, main traits is Ziggurat - most efficient traits to generate random normal samples.</del>
//!     * <del>Code is based on a [C implementation](https://www.seehuhn.de/pages/ziggurat.html) by Jochen Voss.</del>
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         // Normal(mean, std)
//!         let a = Normal(0, 1); // Standard normal
//!         a.sample(100).print();
//!         a.pdf(0).print(); // Maximum probability
//!         a.mean().print();
//!         a.var().print();
//!         a.sd().print();
//!     }
//!     ```
//!
//! ### Beta Distribution
//!
//! ### Gamma Distribution

extern crate rand;
extern crate rand_distr;
use self::rand::distributions::uniform::SampleUniform;
use self::rand::prelude::*;
use self::rand_distr::Distribution;
pub use self::OPDist::*;
pub use self::TPDist::*;
use crate::special::function::*;
//use statistics::rand::ziggurat;
use crate::statistics::stat::Statistics;
use std::convert::Into;
use std::f64::consts::E;

/// One parameter distribution
///
/// # Distributions
/// * `Bernoulli(prob)`: Bernoulli distribution
#[derive(Debug, Clone)]
pub enum OPDist<T: PartialOrd + SampleUniform + Copy + Into<f64>> {
    Bernoulli(T),
    StudentT(T),
}

/// Two parameter distribution
///
/// # Distributions
/// * `Uniform(start, end)`: Uniform distribution
/// * `Normal(mean, std)`: Normal distribution
#[derive(Debug, Clone)]
pub enum TPDist<T: PartialOrd + SampleUniform + Copy + Into<f64>> {
    Uniform(T, T),
    Normal(T, T),
    Beta(T, T),
    Gamma(T, T),
}

/// Extract parameter
pub trait ParametricDist {
    type Parameter;
    fn params(&self) -> Self::Parameter;
}

impl<T: PartialOrd + SampleUniform + Copy + Into<f64>> ParametricDist for OPDist<T> {
    type Parameter = f64;

    fn params(&self) -> Self::Parameter {
        match self {
            Bernoulli(mu) => (*mu).into(),
            StudentT(nu) => (*nu).into(),
        }
    }
}

impl<T: PartialOrd + SampleUniform + Copy + Into<f64>> ParametricDist for TPDist<T> {
    type Parameter = (f64, f64);

    fn params(&self) -> Self::Parameter {
        match self {
            Uniform(a, b) => ((*a).into(), (*b).into()),
            Normal(mu, sigma) => ((*mu).into(), (*sigma).into()),
            Beta(a, b) => ((*a).into(), (*b).into()),
            Gamma(a, b) => ((*a).into(), (*b).into()),
        }
    }
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

    /// Cumulative Distribution Function
    ///
    /// # Type
    /// `f64` -> `f64`
    fn cdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64;
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
            StudentT(nu) => {
                let mut rng = thread_rng();
                let stud = rand_distr::StudentT::<f64>::new((*nu).into()).unwrap();
                stud.sample_iter(&mut rng).take(n).collect()
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
            StudentT(nu) => {
                let dof = (*nu).into();
                let t = x.into();
                1f64 / (dof.sqrt() * beta(0.5f64, dof / 2f64)) * (1f64 + t.powi(2) / dof).powf(-(dof + 1f64) / 2f64)
            }
        }
    }

    fn cdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64 {
        match self {
            Bernoulli(prob) => {
                let k: f64 = x.into();
                if k < 0f64 {
                    0f64
                } else if k < 1f64 {
                    1f64 - (*prob).into()
                } else {
                    1f64
                }
            }
            StudentT(nu) => {
                let x: f64 = x.into();
                let nu: f64 = (*nu).into();
                let _odd_nu = (nu + 1f64) / 2f64;
                let even_nu = nu / 2f64;

                if x > 0f64 {
                    let x_t = nu / (x.powi(2) + nu);
                    1f64 - 0.5 * inc_beta(even_nu, 0.5, x_t)
                } else if x < 0f64 {
                    self.cdf(-x) - 0.5
                } else {
                    0.5
                }
                // 0.5f64 + x * gamma(odd_nu) * hyp2f1(0.5, odd_nu, 1.5, -x.powi(2) / (*nu).into()) / (PI * (*nu).into()).sqrt() * gamma(even_nu)
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
                let normal = rand_distr::Normal::<f64>::new((*m).into(), (*s).into()).unwrap();
                normal.sample_iter(&mut rng).take(n).collect()
            }
//            Normal(m, s) => {
//                let mut rng = thread_rng();
//                let mut v = vec![0f64; n];
//
//                for i in 0..n {
//                    v[i] = ziggurat(&mut rng, (*s).into()) + (*m).into();
//                }
//                v
//            }
            Beta(a, b) => {
                let mut rng = thread_rng();
                let beta = rand_distr::Beta::<f64>::new((*a).into(), (*b).into()).unwrap();
                beta.sample_iter(&mut rng).take(n).collect()
            }
//            Beta(a, b) => {
//                let mut rng1 = thread_rng();
//                let mut rng2 = thread_rng();
//                let mut v = vec![0f64; n];
//
//                let a_f64 = (*a).into();
//                let b_f64 = (*b).into();
//
//                // For acceptance-rejection method
//                let c_x = (a_f64 - 1f64) / (a_f64 + b_f64 - 2f64);
//                let c = self.pdf(c_x); // Beta(mode(x) | a, b)
//
//                let mut iter_num = 0usize;
//
//                while iter_num < n {
//                    let u1 = rng1.gen_range(0f64, 1f64);
//                    let u2 = rng2.gen_range(0f64, 1f64);
//
//                    if u2 <= 1f64 / c * self.pdf(u1) {
//                        v[iter_num] = u1;
//                        iter_num += 1;
//                    }
//                }
//                v
//            }
            Gamma(shape, scale) => {
                let mut rng = thread_rng();
                let gamma = rand_distr::Gamma::<f64>::new((*shape).into(), (*scale).into()).unwrap();
                gamma.sample_iter(&mut rng).take(n).collect()
            }
//            Gamma(a, b) => {
//                let a_f64 = (*a).into();
//                let b_f64 = (*b).into();
//
//                // for Marsaglia & Tsang's Method
//                let d = a_f64 - 1f64 / 3f64;
//                let c = 1f64 / (9f64 * d).sqrt();
//
//                let mut rng1 = thread_rng();
//                let mut rng2 = thread_rng();
//
//                let mut v = vec![0f64; n];
//                let mut iter_num = 0usize;
//
//                while iter_num < n {
//                    let u = rng1.gen_range(0f64, 1f64);
//                    let z = ziggurat(&mut rng2, 1f64);
//                    let w = (1f64 + c * z).powi(3);
//
//                    if z >= -1f64 / c && u.ln() < 0.5 * z.powi(2) + d - d * w + d * w.ln() {
//                        v[iter_num] = d * w / b_f64;
//                        iter_num += 1;
//                    }
//                }
//                v
//            }
        }
    }

    fn pdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64 {
        match self {
            Uniform(a, b) => {
                let val = x.into();
                let a_f64 = (*a).into();
                let b_f64 = (*b).into();
                if val >= a_f64 && val <= b_f64 {
                    let length = b_f64 - a_f64;
                    1f64 / length
                } else {
                    0f64
                }
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
            Gamma(a, b) => {
                let a_f64 = (*a).into();
                let b_f64 = (*b).into();
                1f64 / gamma(a_f64)
                    * b_f64.powf(a_f64)
                    * x.into().powf(a_f64 - 1f64)
                    * E.powf(-b_f64 * x.into())
            }
        }
    }

    fn cdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64 {
        let x: f64 = x.into();
        match self {
            Uniform(a, b) => {
                let a: f64 = (*a).into();
                let b: f64 = (*b).into();

                if x < a {
                    0f64
                } else if x <= b {
                    (x - a) / (b - a)
                } else {
                    1f64
                }
            }
            Normal(m, s) => {
                phi((x - (*m).into()) / (*s).into())
            }
            Beta(a, b) => {
                let a: f64 = (*a).into();
                let b: f64 = (*b).into();

                inc_beta(a, b, x)
            }
            Gamma(a, b) => {
                let a: f64 = (*a).into();
                let b: f64 = (*b).into();

                inc_gamma(a, b * x)
            }
        }
    }
}

impl<T: PartialOrd + SampleUniform + Copy + Into<f64>> Statistics for OPDist<T> {
    type Array = Vec<f64>;
    type Value = f64;

    fn mean(&self) -> Self::Value {
        match self {
            Bernoulli(mu) => (*mu).into(),
            StudentT(_) => 0f64,
        }
    }

    fn var(&self) -> Self::Value {
        match self {
            Bernoulli(mu) => {
                let mu_f64 = (*mu).into();
                mu_f64 * (1f64 - mu_f64)
            }
            StudentT(nu) => {
                let nu_f64 = (*nu).into();
                nu_f64 / (nu_f64 - 2f64)
            }
        }
    }

    fn sd(&self) -> Self::Value {
        match self {
            Bernoulli(_mu) => self.var().sqrt(),
            StudentT(_nu) => self.var().sqrt(),
        }
    }

    fn cov(&self) -> Self::Array {
        unimplemented!()
    }

    fn cor(&self) -> Self::Array {
        unimplemented!()
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
            Gamma(a, b) => (*a).into() / (*b).into(),
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
            Gamma(a, b) => (*a).into() / (*b).into().powi(2),
        }
    }

    fn sd(&self) -> Self::Value {
        match self {
            Uniform(_a, _b) => self.var().sqrt(),
            Normal(_m, s) => (*s).into(),
            Beta(_a, _b) => self.var().sqrt(),
            Gamma(_a, _b) => self.var().sqrt(),
        }
    }

    fn cov(&self) -> Self::Array {
        unimplemented!()
    }

    fn cor(&self) -> Self::Array {
        unimplemented!()
    }
}
