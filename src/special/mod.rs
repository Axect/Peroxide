//! Special functions module
//!
//! This module provides implementations of various special mathematical functions
//! commonly used in statistical and scientific computing. It includes:
//!
//! - Basic special functions:
//!   - Gaussian (Normal) function
//!   - Gamma function and its logarithm
//!   - Pochhammer symbol (rising factorial)
//!
//! - Incomplete special functions:
//!   - Regularized incomplete gamma function and its inverse
//!   - Regularized incomplete beta function and its inverse
//!
//! - Error functions:
//!   - Error function (erf) and its complement (erfc)
//!   - Inverse error function and inverse complementary error function
//!
//! - Other functions:
//!   - Beta function
//!   - Phi function (CDF of the standard normal distribution)
//!   - Lambert W function (principal branch W₀ and secondary branch W₋₁)
//!
//! Many of these functions are implemented using efficient numerical approximations
//! or by wrapping functions from other crates (e.g., `puruspe`).
//!
//! The module also includes an enum `LambertWAccuracyMode` to control the
//! accuracy-speed trade-off for Lambert W function calculations.

pub mod function;
pub mod lanczos;
