extern crate special;
extern crate special_fun;
use self::special::{Error, Gamma};

use std::f64::consts::PI;
use special::lanczos::{gamma_approx, ln_gamma_approx};

/// Gaussian function
///
/// `N(x|μ,σ) = 1/√(2πσ^2) exp(-(x-μ)^2/(2σ^2))`
pub fn gaussian(x: f64, mu: f64, sigma: f64) -> f64 {
    1f64 / ((2f64 * PI).sqrt() * sigma) * (-0.5 * ((x - mu) / sigma).powi(2)).exp()
}

/// Gamma function
///
/// # Description
/// Use Lanczos approximation to implement Gamma function ($g=5, n=7$)
///
/// # References
/// * [Robert Munafo, Coefficients for the Lanczos Approximation to the Gamma Function](https://mrob.com/pub/ries/lanczos-gamma.html)
/// * [Paul Godfrey, A note on the computation of the convergent Lanczos complex Gamma approximation (web page), 2001.](http://my.fit.edu/~gabdo/gamma.txt)
pub fn gamma(x: f64) -> f64 {
    gamma_approx(x)
}

/// Logarithm Gamma function
///
/// # Description
/// Use Lanczos approximation to implement Gamma function ($g=5, n=7$)
///
/// # References
/// * [Robert Munafo, Coefficients for the Lanczos Approximation to the Gamma Function](https://mrob.com/pub/ries/lanczos-gamma.html)
/// * [Paul Godfrey, A note on the computation of the convergent Lanczos complex Gamma approximation (web page), 2001.](http://my.fit.edu/~gabdo/gamma.txt)
pub fn ln_gamma(x: f64) -> f64 {
    ln_gamma_approx(x)
}

/// Pochhammer symbol
pub fn poch(x: f64, n: usize) -> f64 {
    let mut s = 1f64;
    for i in 0 .. n {
        s *= x + i as f64;
    }
    s
}

/// Digamma function
///
/// Wrapper of `digamma` function of `special` crate
pub fn digamma(x: f64) -> f64 {
    x.digamma()
}

/// Regularized incomplete gamma integral (Lower)
///
/// Wrapper of `gammp` function of `puruspe` crate
pub fn inc_gamma(a: f64, x: f64) -> f64 {
    puruspe::gammp(a, x)
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