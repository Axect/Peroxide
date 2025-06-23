//! Probabilistic distributions
//!
//! ## Probability Distribution
//!
//! * There are some famous pdf in Peroxide (not checked pdfs will be implemented soon)
//!     * Bernoulli
//!     * Binomial
//!     * Beta
//!     * Dirichlet
//!     * Gamma
//!     * Normal
//!     * Student's t
//!     * Uniform
//!     * Weighted Uniform
//!     * Log Normal
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
//!     * `sample_with_rng`: Extract samples with specific rng
//!     * `pdf` : Calculate pdf value at specific point
//!     ```no_run
//!     use rand::Rng;
//!     use rand_distr::uniform::SampleUniform;
//!     pub trait RNG {
//!         /// Extract samples of distributions
//!         fn sample(&self, n: usize) -> Vec<f64>;
//!
//!         /// Extract samples of distributions with rng
//!         fn sample_with_rng<R: Rng>(&self, rng: &mut R, n: usize) -> Vec<f64>;
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
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         let mut rng = smallrng_from_seed(42);
//!         let b = Bernoulli(0.1);                   // Bern(x | 0.1)
//!         b.sample(100).print();                    // Generate 100 samples
//!         b.sample_with_rng(&mut rng, 100).print(); // Generate 100 samples with specific rng
//!         b.pdf(0).print();                         // 0.9
//!         b.mean().print();                         // 0.1
//!         b.var().print();                          // 0.09 (approximately)
//!         b.sd().print();                           // 0.3  (approximately)
//!     }
//!     ```
//!
//! ### Uniform Distribution
//!
//! * Definition
//!     $$\text{Unif}(x | a, b) = \begin{cases} \frac{1}{b - a} & x \in \[a,b\]\\\ 0 & \text{otherwise} \end{cases}$$
//! * Representative value
//!     * Mean: $\frac{a + b}{2}$
//!     * Var : $\frac{1}{12}(b-a)^2$
//! * To generate uniform random number, Peroxide uses `rand` crate
//! * **Caution**: `Uniform(T, T)` generates `T` type samples (only for `Uniform`)
//!
//!     ```rust
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         // Uniform(start, end)
//!         let a = Uniform(0f64, 1f64); // It will generate `f64` samples.
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

//! ### Beta Distribution
//!
//! * Definition
//!     $$\text{Beta}(x | \alpha, \beta) = \frac{1}{\text{B}(\alpha, \beta)} x^{\alpha-1} (1-x)^{\beta-1}$$
//!     where $\text{B}(\alpha, \beta) = \frac{\Gamma(\alpha)\Gamma(\beta)}{\Gamma(\alpha+\beta)}$ is the Beta function.
//! * Representative value
//!     * Mean: $\frac{\alpha}{\alpha+\beta}$
//!     * Var: $\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}$
//! * To generate beta random samples, Peroxide uses the `rand_distr::Beta` distribution from the `rand_distr` crate.
//!
//!     ```rust
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         // Beta(alpha, beta)
//!         let a = Beta(2.0, 5.0);
//!         a.sample(100).print();
//!         a.pdf(0.3).print();
//!         a.mean().print();
//!         a.var().print();
//!     }
//!     ```
//!
//! ### Gamma Distribution
//!
//! * Definition
//!     $$\text{Gamma}(x | \alpha, \beta) = \frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha-1} e^{-\beta x}$$
//!     where $\Gamma(\alpha) = \int_0^\infty x^{\alpha-1} e^{-x} dx$ is the Gamma function.
//! * Representative value
//!     * Mean: $\frac{\alpha}{\beta}$
//!     * Var: $\frac{\alpha}{\beta^2}$
//! * To generate gamma random samples, Peroxide uses the `rand_distr::Gamma` distribution from the `rand_distr` crate.
//!
//!     ```rust
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         // Gamma(shape, scale)
//!         let a = Gamma(2.0, 1.0);
//!         a.sample(100).print();
//!         a.pdf(1.5).print();
//!         a.mean().print();
//!         a.var().print();
//!     }
//!     ```
//!
//! ### Binomial Distribution
//!
//! * Definition
//!     $$\text{Binom}(k | n, p) = \binom{n}{k} p^k (1-p)^{n-k}$$
//!     where $\binom{n}{k} = \frac{n!}{k!(n-k)!}$ is the binomial coefficient.
//! * Representative value
//!     * Mean: $np$
//!     * Var: $np(1-p)$
//! * To generate binomial random samples, Peroxide uses the `rand_distr::Binomial` distribution from the `rand_distr` crate.
//!
//!     ```rust
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         // Binomial(n, p)
//!         let a = Binomial(10, 0.3);
//!         a.sample(100).print();
//!         a.pdf(3).print();
//!         a.mean().print();
//!         a.var().print();
//!     }
//!     ```
//!
//! ### Student's t Distribution
//!
//! * Definition
//!     $$\text{StudentT}(x | \nu) = \frac{\Gamma(\frac{\nu+1}{2})}{\sqrt{\nu\pi}\,\Gamma(\frac{\nu}{2})} \left(1+\frac{x^2}{\nu} \right)^{-\frac{\nu+1}{2}}$$
//!     where $\nu$ is the degrees of freedom and $\Gamma$ is the Gamma function.
//! * Representative value
//!     * Mean: 0 (for $\nu > 1$)
//!     * Var: $\frac{\nu}{\nu-2}$ (for $\nu > 2$)
//! * To generate Student's t random samples, Peroxide uses the `rand_distr::StudentT` distribution from the `rand_distr` crate.
//!
//!     ```rust
//!     use peroxide::fuga::*;
//!
//!     fn main() {
//!         // StudentT(nu)
//!         let a = StudentT(5.0);
//!         a.sample(100).print();
//!         a.pdf(1.0).print();
//!         a.mean().print(); // Undefined for nu <= 1
//!         a.var().print();  // Undefined for nu <= 2
//!     }
//!     ```
//!
//! ### Weighted Uniform Distribution
//!
//! * Definition
//!    $$\text{WUnif}(x | \mathbf{W}, \mathcal{I}) = \frac{1}{\sum_{j=1}^n w_j \mu(I_j)} \sum_{i=1}^n w_i
//!    \mathbb{1}_{I_i}(x)$$
//!    * $\mathbf{W} = (w_i)$: Weights
//!    * $\mathcal{I} = \\{I_i\\}$: Intervals
//!    * $\mu(I_i)$: Measure of $I_i$
//!    * $\mathbb{1}_{I_i}(x)$: Indicator function
//!
//! * Reference
//!     * [Piecewise Rejection Sampling](https://axect.github.io/posts/006_prs/#22-weighted-uniform-distribution)
//!
//! ### Log Normal Distribution
//!
//! * Definition
//!     $$\text{LogNormal}(x | \mu, \sigma) = \frac{1}{x\sigma\sqrt{2\pi}} e^{-\frac{(\ln x -
//!     \mu)^2}{2\sigma^2}}$$
//!     where $\mu$ is the mean of the natural logarithm of the variable and $\sigma$ is the
//!     standard deviation of the natural logarithm of the variable.
//! * Representative value
//!     * Mean: $e^{\mu + \frac{\sigma^2}{2}}$
//!     * Var: $(e^{\sigma^2} - 1)e^{2\mu + \sigma^2}$
//! * To generate log-normal random samples, Peroxide uses the `rand_distr::LogNormal` distribution from the `rand_distr` crate.

extern crate rand;
extern crate rand_distr;
use rand_distr::weighted::WeightedAliasIndex;

use self::rand::prelude::*;
use self::rand_distr::uniform::SampleUniform;
pub use self::OPDist::*;
pub use self::TPDist::*;
use crate::special::function::*;
use crate::traits::fp::FPVector;
//use statistics::rand::ziggurat;
use self::WeightedUniformError::*;
use crate::statistics::{ops::C, stat::Statistics};
use crate::util::non_macro::{linspace, seq};
use crate::util::useful::{auto_zip, find_interval};
use anyhow::{bail, Result};
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
    Binomial(usize, T),
    Normal(T, T),
    Beta(T, T),
    Gamma(T, T),
    LogNormal(T, T),
}

pub struct WeightedUniform<T: PartialOrd + SampleUniform + Copy + Into<f64>> {
    weights: Vec<T>,
    sum: T,
    intervals: Vec<(T, T)>,
}

#[derive(Debug, Clone, Copy)]
pub enum WeightedUniformError {
    AllZeroWeightError,
    LengthMismatchError,
    NoNonZeroIntervalError,
    EmptyWeightError,
}

impl std::fmt::Display for WeightedUniformError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WeightedUniformError::AllZeroWeightError => write!(f, "all weights are zero"),
            WeightedUniformError::LengthMismatchError => {
                write!(f, "weights and intervals have different length")
            }
            WeightedUniformError::NoNonZeroIntervalError => write!(f, "no non-zero interval found"),
            WeightedUniformError::EmptyWeightError => write!(f, "weights are empty"),
        }
    }
}

impl WeightedUniform<f64> {
    /// Create a new weighted uniform distribution
    ///
    /// # Examples
    /// ```
    /// use peroxide::fuga::*;
    ///
    /// fn main() -> Result<(), Box<dyn Error>> {
    ///     let weights = vec![1f64, 3f64, 0f64, 2f64];
    ///     let intervals = vec![0f64, 1f64, 2f64, 4f64, 5f64];
    ///     let w = WeightedUniform::new(weights, intervals)?;
    ///     assert_eq!(w.weights(), &vec![1f64, 3f64, 2f64]);
    ///     assert_eq!(w.intervals(), &vec![(0f64, 1f64), (1f64, 2f64), (4f64, 5f64)]);
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn new(weights: Vec<f64>, intervals: Vec<f64>) -> Result<Self> {
        let mut weights = weights;
        if weights.len() == 0 {
            bail!(EmptyWeightError);
        }
        if weights.iter().all(|&x| x == 0f64) {
            bail!(AllZeroWeightError);
        }
        let mut intervals = auto_zip(&intervals);
        if weights.len() != intervals.len() {
            bail!(LengthMismatchError);
        }

        // Remove zero weights & corresponding intervals
        let mut i = 0;
        let mut j = weights.len();
        while i < j {
            if weights[i] == 0f64 {
                weights.remove(i);
                intervals.remove(i);
                j -= 1;
            } else {
                i += 1;
            }
        }

        let sum = weights
            .iter()
            .zip(intervals.iter())
            .fold(0f64, |acc, (w, (a, b))| acc + w * (b - a));

        Ok(WeightedUniform {
            weights,
            sum,
            intervals,
        })
    }

    /// Create WeightedUniform from max pooling
    ///
    /// # Examples
    /// ```
    /// use peroxide::fuga::*;
    ///
    /// fn main() -> Result<(), Box<dyn Error>> {
    ///     let w = WeightedUniform::from_max_pool_1d(f, (-2f64, 3f64), 10, 1e-3)?;
    ///     w.weights().print();
    ///
    ///     Ok(())
    /// }
    ///
    /// fn f(x: f64) -> f64 {
    ///     if x.abs() < 1f64 {
    ///         1f64 - x.abs()
    ///     } else {
    ///        0f64
    ///     }
    /// }
    /// ```
    pub fn from_max_pool_1d<F>(f: F, (a, b): (f64, f64), n: usize, eps: f64) -> Result<Self>
    where
        F: Fn(f64) -> f64 + Copy,
    {
        // Find non-zero intervals
        let mut a = a;
        let mut b = b;
        let trial = seq(a, b, eps);
        for i in 0..trial.len() {
            let x = trial[i];
            if f(x) > 0f64 {
                a = if i > 0 { trial[i - 1] } else { x };
                break;
            }
        }
        for i in (0..trial.len()).rev() {
            let x = trial[i];
            if f(x) > 0f64 {
                b = if i < trial.len() - 1 { trial[i + 1] } else { x };
                break;
            }
        }
        if a >= b {
            bail!(NoNonZeroIntervalError);
        }
        let domain = linspace(a, b, n + 1);

        // Find intervals
        let intervals = auto_zip(&domain);

        // Find weights
        let weights: Vec<f64> = intervals
            .iter()
            .map(|(a, b)| seq(*a, *b + eps, eps).reduce(0f64, |acc, x| acc.max(f(x))))
            .collect();

        Self::new(weights, domain)
    }

    pub fn weights(&self) -> &Vec<f64> {
        &self.weights
    }

    pub fn intervals(&self) -> &Vec<(f64, f64)> {
        &self.intervals
    }

    pub fn domain_linspace(&self, n: usize) -> Vec<f64> {
        linspace(
            self.intervals[0].0,
            self.intervals[self.intervals.len() - 1].1,
            n,
        )
    }

    pub fn domain_seq(&self, step: f64) -> Vec<f64> {
        seq(
            self.intervals[0].0,
            self.intervals[self.intervals.len() - 1].1,
            step,
        )
    }

    pub fn sum(&self) -> f64 {
        self.sum
    }

    pub fn update_weights(&mut self, weights: Vec<f64>) {
        assert_eq!(self.intervals.len(), weights.len());
        self.weights = weights;
        self.sum = self
            .weights
            .iter()
            .zip(self.intervals.iter())
            .fold(0f64, |acc, (w, (a, b))| acc + w * (b - a));
    }

    pub fn update_intervals(&mut self, intervals: Vec<f64>) {
        assert_eq!(self.weights.len() + 1, intervals.len());
        self.intervals = auto_zip(&intervals);
        self.sum = self
            .weights
            .iter()
            .zip(self.intervals.iter())
            .fold(0f64, |acc, (w, (a, b))| acc + w * (b - a));
    }

    pub fn weight_at(&self, x: f64) -> f64 {
        let i = find_interval(self.intervals(), x);
        self.weights[i]
    }

    pub fn interval_at(&self, x: f64) -> (f64, f64) {
        let i = find_interval(self.intervals(), x);
        self.intervals[i]
    }
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
            Binomial(a, b) => (*a as f64, (*b).into()),
            Normal(mu, sigma) => ((*mu).into(), (*sigma).into()),
            Beta(a, b) => ((*a).into(), (*b).into()),
            Gamma(a, b) => ((*a).into(), (*b).into()),
            LogNormal(mu, sigma) => ((*mu).into(), (*sigma).into()),
        }
    }
}

impl<T: PartialOrd + SampleUniform + Copy + Into<f64>> ParametricDist for WeightedUniform<T> {
    type Parameter = (Vec<f64>, Vec<(f64, f64)>);

    fn params(&self) -> Self::Parameter {
        let weights = self.weights.iter().map(|x| (*x).into()).collect();
        let intervals = self
            .intervals
            .iter()
            .map(|(x, y)| ((*x).into(), (*y).into()))
            .collect();
        (weights, intervals)
    }
}

/// Random Number Generator trait
///
/// # Methods
/// * `sample`: extract samples
pub trait RNG {
    /// Extract samples of distributions
    fn sample(&self, n: usize) -> Vec<f64> {
        let mut rng = rand::rng();
        self.sample_with_rng(&mut rng, n)
    }

    /// Extract samples of distributions with rng
    fn sample_with_rng<R: Rng + Clone>(&self, rng: &mut R, n: usize) -> Vec<f64>;

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
    fn sample_with_rng<R: Rng + Clone>(&self, rng: &mut R, n: usize) -> Vec<f64> {
        match self {
            Bernoulli(prob) => {
                assert!(
                    (*prob).into() <= 1f64,
                    "Probability should be smaller than 1"
                );

                let mut v = vec![0f64; n];

                for i in 0..n {
                    let uniform = rng.random_range(0f64..=1f64);
                    if uniform <= (*prob).into() {
                        v[i] = 1f64;
                    } else {
                        v[i] = 0f64;
                    }
                }
                v
            }
            StudentT(nu) => {
                let stud = rand_distr::StudentT::<f64>::new((*nu).into()).unwrap();
                stud.sample_iter(rng).take(n).collect()
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
                1f64 / (dof.sqrt() * beta(0.5f64, dof / 2f64))
                    * (1f64 + t.powi(2) / dof).powf(-(dof + 1f64) / 2f64)
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
    fn sample_with_rng<R: Rng + Clone>(&self, rng: &mut R, n: usize) -> Vec<f64> {
        match self {
            Uniform(start, end) => {
                let mut v = vec![0f64; n];

                for i in 0..n {
                    v[i] = rng.random_range(*start..=*end).into();
                }
                v
            }

            Binomial(num, mu) => {
                let binom = rand_distr::Binomial::new(*num as u64, (*mu).into()).unwrap();
                binom.sample_iter(rng).take(n).map(|t| t as f64).collect()
            }

            Normal(m, s) => {
                let normal = rand_distr::Normal::<f64>::new((*m).into(), (*s).into()).unwrap();
                normal.sample_iter(rng).take(n).collect()
            }
            //            Normal(m, s) => {
            //                let mut rng = rand::rng();
            //                let mut v = vec![0f64; n];
            //
            //                for i in 0..n {
            //                    v[i] = ziggurat(&mut rng, (*s).into()) + (*m).into();
            //                }
            //                v
            //            }
            Beta(a, b) => {
                let beta = rand_distr::Beta::<f64>::new((*a).into(), (*b).into()).unwrap();
                beta.sample_iter(rng).take(n).collect()
            }
            //            Beta(a, b) => {
            //                let mut rng1 = rand::rng();
            //                let mut rng2 = rand::rng();
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
            //                    let u1 = rng1.random_range(0f64, 1f64);
            //                    let u2 = rng2.random_range(0f64, 1f64);
            //
            //                    if u2 <= 1f64 / c * self.pdf(u1) {
            //                        v[iter_num] = u1;
            //                        iter_num += 1;
            //                    }
            //                }
            //                v
            //            }
            Gamma(shape, scale) => {
                let gamma =
                    rand_distr::Gamma::<f64>::new((*shape).into(), (*scale).into()).unwrap();
                gamma.sample_iter(rng).take(n).collect()
            } //            Gamma(a, b) => {
              //                let a_f64 = (*a).into();
              //                let b_f64 = (*b).into();
              //
              //                // for Marsaglia & Tsang's Method
              //                let d = a_f64 - 1f64 / 3f64;
              //                let c = 1f64 / (9f64 * d).sqrt();
              //
              //                let mut rng1 = rand::rng();
              //                let mut rng2 = rand::rng();
              //
              //                let mut v = vec![0f64; n];
              //                let mut iter_num = 0usize;
              //
              //                while iter_num < n {
              //                    let u = rng1.random_range(0f64, 1f64);
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
            LogNormal(mu, sigma) => {
                let log_normal =
                    rand_distr::LogNormal::<f64>::new((*mu).into(), (*sigma).into()).unwrap();
                log_normal.sample_iter(rng).take(n).collect()
            }
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
            Binomial(n, mu) => {
                let n = *n;
                let mu = (*mu).into();
                let m = x.into() as usize;
                (C(n, m) as f64) * mu.powi(m as i32) * (1f64 - mu).powi((n - m) as i32)
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
            LogNormal(mu, sigma) => {
                let mu_f64 = (*mu).into();
                let sigma_f64 = (*sigma).into();
                if x.into() <= 0f64 {
                    0f64
                } else {
                    1f64 / (x.into() * sigma_f64 * (2f64 * std::f64::consts::PI).sqrt())
                        * E.powf(-((x.into().ln() - mu_f64).powi(2)) / (2f64 * sigma_f64.powi(2)))
                }
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
            Binomial(n, mu) => {
                let n = *n;
                let p = (*mu).into();
                let q = 1f64 - p;
                let k: f64 = x.into();
                inc_beta(n as f64 - k, k + 1f64, q)
            }
            Normal(m, s) => phi((x - (*m).into()) / (*s).into()),
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
            LogNormal(mu, sigma) => {
                phi(
                    (x.ln() - (*mu).into()) / (*sigma).into(),
                )
            }
        }
    }
}

impl RNG for WeightedUniform<f64> {
    fn sample_with_rng<R: Rng + Clone>(&self, rng: &mut R, n: usize) -> Vec<f64> {
        let w = WeightedAliasIndex::new(self.weights.clone()).unwrap();
        let mut rng_clip = rng.clone();
        let ics: Vec<usize> = w.sample_iter(&mut rng_clip).take(n).collect();
        *rng = rng_clip;

        ics.into_iter()
            .map(|idx| {
                let (l, r) = self.intervals[idx];
                rng.random_range(l..=r)
            })
            .collect::<Vec<f64>>()
    }

    fn pdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64 {
        let x: f64 = x.into();
        if x < self.intervals[0].0 || x > self.intervals[self.intervals.len() - 1].1 {
            return 0f64;
        }
        let idx = find_interval(self.intervals(), x);
        self.weights[idx] / self.sum
    }

    fn cdf<S: PartialOrd + SampleUniform + Copy + Into<f64>>(&self, x: S) -> f64 {
        let x: f64 = x.into();
        if x < self.intervals[0].0 {
            return 0f64;
        } else if x > self.intervals[self.intervals.len() - 1].1 {
            return 1f64;
        }
        let idx = find_interval(self.intervals(), x);
        self.weights[0..=idx]
            .iter()
            .zip(self.intervals[0..=idx].iter())
            .fold(0f64, |acc, (w, (a, b))| acc + w * (b - a))
            / self.sum
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
            Binomial(n, mu) => (*n as f64) * (*mu).into(),
            Normal(m, _s) => (*m).into(),
            Beta(a, b) => (*a).into() / ((*a).into() + (*b).into()),
            Gamma(a, b) => (*a).into() / (*b).into(),
            LogNormal(mu, sigma) => {
                let mu_f64 = (*mu).into();
                let sigma_f64 = (*sigma).into();
                E.powf(mu_f64 + 0.5 * sigma_f64.powi(2))
            }
        }
    }

    fn var(&self) -> Self::Value {
        match self {
            Uniform(a, b) => ((*b).into() - (*a).into()).powi(2) / 12f64,
            Binomial(n, mu) => (*n as f64) * (*mu).into() * (1f64 - (*mu).into()),
            Normal(_m, s) => (*s).into().powi(2),
            Beta(a, b) => {
                let a_f64 = (*a).into();
                let b_f64 = (*b).into();
                a_f64 * b_f64 / ((a_f64 + b_f64).powi(2) * (a_f64 + b_f64 + 1f64))
            }
            Gamma(a, b) => (*a).into() / (*b).into().powi(2),
            LogNormal(mu, sigma) => {
                let mu_f64 = (*mu).into();
                let sigma_f64 = (*sigma).into();
                (E.powf(sigma_f64.powi(2)) - 1f64) * E.powf(2f64 * mu_f64 + sigma_f64.powi(2))
            }
        }
    }

    fn sd(&self) -> Self::Value {
        match self {
            Uniform(_a, _b) => self.var().sqrt(),
            Binomial(_n, _mu) => self.var().sqrt(),
            Normal(_m, s) => (*s).into(),
            Beta(_a, _b) => self.var().sqrt(),
            Gamma(_a, _b) => self.var().sqrt(),
            LogNormal(_mu, _sigma) => self.var().sqrt(),
        }
    }

    fn cov(&self) -> Self::Array {
        unimplemented!()
    }

    fn cor(&self) -> Self::Array {
        unimplemented!()
    }
}

impl Statistics for WeightedUniform<f64> {
    type Array = Vec<f64>;
    type Value = f64;

    fn mean(&self) -> Self::Value {
        self.intervals()
            .iter()
            .zip(self.weights().iter())
            .map(|((l, r), w)| (r.powi(2) - l.powi(2)) / 2f64 * w)
            .sum::<f64>()
            / self.sum
    }

    fn var(&self) -> Self::Value {
        let mean = self.mean();
        self.intervals()
            .iter()
            .zip(self.weights().iter())
            .map(|((l, r), w)| w * (r.powi(3) - l.powi(3)) / 3f64)
            .sum::<f64>()
            / self.sum
            - mean * mean
    }

    fn sd(&self) -> Self::Value {
        self.var().sqrt()
    }

    fn cov(&self) -> Self::Array {
        vec![self.var()]
    }

    fn cor(&self) -> Self::Array {
        vec![1f64]
    }
}
