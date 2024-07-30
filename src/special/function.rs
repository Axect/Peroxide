use crate::special::lanczos::{gamma_approx, ln_gamma_approx};
use std::f64::consts::PI;

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
    for i in 0..n {
        s *= x + i as f64;
    }
    s
}

// /// Digamma function
// ///
// /// Wrapper of `digamma` function of `special` crate
// pub fn digamma(x: f64) -> f64 {
//     x.digamma()
// }

/// Regularized incomplete gamma integral (Lower)
///
/// Wrapper of `gammp` function of `puruspe` crate
pub fn inc_gamma(a: f64, x: f64) -> f64 {
    puruspe::gammp(a, x)
}

/// Inverse of regularized incomplete gamma integral (Lower)
///
/// Wrapper of `invgammp` function of `puruspe` crate
pub fn inv_inc_gamma(p: f64, a: f64) -> f64 {
    puruspe::invgammp(p, a)
}

/// Error function
///
/// Wrapper of `erf` function of `puruspe` crate
pub fn erf(x: f64) -> f64 {
    puruspe::erf(x)
}

/// Complement error function
///
/// Wrapper of `erfc` function of `puruspe` crate
pub fn erfc(x: f64) -> f64 {
    puruspe::erfc(x)
}

/// Inverse error function
///
/// Wrapper of `inverf` function of `puruspe` crate
pub fn inv_erf(x: f64) -> f64 {
    puruspe::inverf(x)
}

/// Inverse complementary error function
///
/// Wrapper of `inverfc` function of `puruspe` crate
pub fn inv_erfc(p: f64) -> f64 {
    puruspe::inverfc(p)
}

/// Beta function
///
/// Wrapper of `beta` function of `puruspe` crate
pub fn beta(a: f64, b: f64) -> f64 {
    puruspe::beta(a, b)
}

/// Regularized incomplete Beta function
///
/// Wrapper of `betai` function of `puruspe` crate
pub fn inc_beta(a: f64, b: f64, x: f64) -> f64 {
    puruspe::betai(a, b, x)
}

/// Inverse regularized incomplete beta function
///
/// Wrapper of `invbetai` function of `puruspe` crate
pub fn inv_inv_beta(p: f64, a: f64, b: f64) -> f64 {
    puruspe::invbetai(p, a, b)
}

/// Phi (CDF for Normal Dist)
///
/// $$\Phi(x) = \frac{1}{2}\left[1 + \text{erf}\left(\frac{x}{\sqrt{2}}\right) \right]$$
pub fn phi(x: f64) -> f64 {
    0.5 * (1f64 + erf(x / 2f64.sqrt()))
}

/// The principal branch of the Lambert W function, W_0(`z`).
/// 
/// Returns [`NAN`](f64::NAN) if the given input is smaller than -1/e.
/// 
/// Wrapper of `lambert_w_0` function of `lambert_w` crate.
pub fn lambert_w0(z: f64) -> f64 {
    lambert_w::lambert_w_0(z).unwrap_or(f64::NAN)
}

/// The secondary branch of the Lambert W function, W_-1(`z`).
/// 
/// Returns [`NAN`](f64::NAN) if the given input is smaller than -1/e or positive.
/// 
/// Wrapper of `lambert_w_m1` function of `lambert_w` crate.
pub fn lambert_wm1(z: f64) -> f64 {
    lambert_w::lambert_w_m1(z).unwrap_or(f64::NAN)
}

// /// Hypergeometric function 2F1
// ///
// /// Wrapper of `hyp2f1` function of `special-fun` crate
// pub fn hyp2f1(a: f64, b: f64, c: f64, x: f64) -> f64 {
//     unsafe {
//         special_fun::unsafe_cephes_double::hyp2f1(a, b, c, x)
//     }
// }
