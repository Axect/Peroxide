use std::f64::consts::{PI, E};

/// Gaussian function
///
/// `N(x|μ,σ) = 1/√(2πσ^2) exp(-(x-μ)^2/(2σ^2))`
pub fn gaussian(x: f64, mu: f64, sigma: f64) -> f64 {
    1f64/((2f64*PI).sqrt()*sigma) * (- 0.5 * ((x - mu)/sigma).powi(2)).exp()
}
