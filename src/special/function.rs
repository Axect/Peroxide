extern crate special;
extern crate special_fun;
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
pub fn beta(a: f64, b: f64) -> f64 {
    special_fun::cephes_double::beta(a, b)
}

/// Incomplete Beta function
///
/// Wrapper of `incbet` function of `special-fun` crate
pub fn inc_beta(a: f64, b: f64, x: f64) -> f64 {
    special_fun::cephes_double::incbet(a, b, x)
}

/// Phi (CDF for Normal Dist)
///
/// $$\Phi(x) = \frac{1}{2}\left[1 + \text{erf}\left(\frac{x}{\sqrt{2}}\right) \right]$$
pub fn phi(x: f64) -> f64 {
    0.5 * (1f64 + erf(x / 2f64.sqrt()))
}

/// Hypergeometric function 2F1
///
/// Wrapper of `hyp2f1` function of `special-fun` crate
pub fn hyp2f1(a: f64, b: f64, c: f64, x: f64) -> f64 {
    unsafe {
        special_fun::unsafe_cephes_double::hyp2f1(a, b, c, x)
    }
}