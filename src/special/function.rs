extern crate special;
use self::special::{Error, Gamma};

use std::f64::consts::PI;

/// Gaussian function
///
/// `N(x|μ,σ) = 1/√(2πσ^2) exp(-(x-μ)^2/(2σ^2))`
pub fn gaussian(x: f64, mu: f64, sigma: f64) -> f64 {
    1f64 / ((2f64 * PI).sqrt() * sigma) * (-0.5 * ((x - mu) / sigma).powi(2)).exp()
}

/// Gamma function
///
/// Wrapper of `gamma` function of `special` crate
pub fn gamma(x: f64) -> f64 {
    x.gamma()
}

/// Digamma function
///
/// Wrapper of `digamma` function of `special` crate
pub fn digamma(x: f64) -> f64 {
    x.digamma()
}

/// Error function
///
/// Wrapper of `error` function of `special` crate
pub fn erf(x: f64) -> f64 {
    x.error()
}

/// Complement error function
///
/// Wrapper of `compl_error` function of `special` crate
pub fn erfc(x: f64) -> f64 {
    x.compl_error()
}

/// Inverse error function
///
/// Wrapper of `inv_error` function of `special` crate
pub fn erf_inv(x: f64) -> f64 {
    x.inv_error()
}

/// Beta function
///
/// Wrapper of `inc_beta` function of `special` crate
pub fn beta(a: f64, b: f64) -> f64 {
    gamma(a) * gamma(b) / gamma(a + b)
}
