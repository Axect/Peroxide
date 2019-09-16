//! `peroxide` is comprehensive numerical library for Rust.
//!
//! ## Components
//!
//! `peroxide` has various components for scientific computation.
//!
//! * Linear Algebra
//!     * Matrix operations
//!         * `+,-,*,/`
//!         * LU, Determinant, Inverse
//!     * Vector operations
//! * Statistics
//!     * Statistical operations
//!         * `mean, var, sd`
//!         * `factorial, P, C, H`
//!     * Distributions
//!         * Bernoulli
//!         * Uniform
//!         * Normal
//!         * Gamma
//!         * Beta
//! * Special functions (Using `special` crate)
//!     * Gaussian
//!     * Gamma
//!     * Beta
//!     * Error
//! * Automatic Differentiation
//!     * Dual number system
//!     * Hyper dual number system
//! * Numerical Utils
//!     * Interpolation
//!     * Spline
//!     * Polynomial
//! * Differential Equations
//!     * Explicit
//!         * Runge-Kutta 4th order
//!         * Euler methods
//!     * Implicit
//!         * Backward Euler
//!         * Gauss-Legendre 4th order
//! * Communication with Python
//!     * Support `pickle` type
//!     * Plot with `matplotlib`
//!
//! ## Quick Start
//!
//! ### Cargo.toml
//!
//! * To use `peroxide`, you should edit `Cargo.toml`
//! * Current document version is corresponding to `0.11.5`
//!
//!     ```toml
//!     peroxide = "0.11"
//!     ```
//!
//! ## Import all at once
//!
//! * You can import all functions & structures at once
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::*;
//!
//!     fn main() {
//!         //Some codes...
//!     }
//!     ```

#[cfg(feature = "oxidize")]
extern crate blas;

#[cfg(feature = "oxidize")]
extern crate lapack;

#[cfg(feature = "plot")]
extern crate pyo3;

#[cfg(feature = "simd")]
extern crate packed_simd;

#[cfg(feature = "serde")]
extern crate serde;

extern crate rand;

pub mod statistics;
pub mod structure;
#[macro_use]
pub mod macros;
pub mod ml;
pub mod numerical;
pub mod operation;
pub mod redox;
pub mod special;
pub mod util;

#[allow(unused_imports)]
pub use macros::{julia_macro::*, matlab_macro::*, r_macro::*};

#[allow(unused_imports)]
pub use structure::matrix::*;

#[allow(unused_imports)]
pub use structure::vector::*;

#[allow(unused_imports)]
pub use statistics::stat::*;

#[allow(unused_imports)]
pub use macros::r_macro::*;

#[allow(unused_imports)]
pub use macros::matlab_macro::*;

#[allow(unused_imports)]
pub use macros::julia_macro::*;

#[allow(unused_imports)]
pub use statistics::rand::*;

#[allow(unused_imports)]
pub use util::print::*;

#[allow(unused_imports)]
pub use util::non_macro::*;

#[allow(unused_imports)]
pub use structure::polynomial::*;

#[allow(unused_imports)]
pub use numerical::interp::*;

#[allow(unused_imports)]
pub use numerical::spline::*;

#[allow(unused_imports)]
pub use ml::reg::*;

#[allow(unused_imports)]
pub use structure::dual::*;

#[allow(unused_imports)]
pub use operation::extra_ops::*;

#[allow(unused_imports)]
pub use util::useful::*;

#[allow(unused_imports)]
pub use structure::multinomial::*;

#[allow(unused_imports)]
pub use numerical::utils::*;

//#[allow(unused_imports)]
//pub use numerical::newton::*;

//#[allow(unused_imports)]
//pub use numerical::bdf::*;

#[allow(unused_imports)]
pub use util::api::*;

//#[allow(unused_imports)]
//pub use numerical::gauss_legendre::*;

#[allow(unused_imports)]
pub use statistics::dist::*;

#[allow(unused_imports)]
pub use special::function::*;

#[allow(unused_imports)]
pub use statistics::ops::*;

#[allow(unused_imports)]
pub use structure::hyper_dual::*;

#[allow(unused_imports)]
pub use util::writer::*;

#[allow(unused_imports)]
pub use operation::mut_ops::*;

#[allow(unused_imports)]
pub use numerical::ode::*;

#[allow(unused_imports)]
pub use operation::number::*;

#[allow(unused_imports)]
#[cfg(feature = "plot")]
pub use util::plot::*;

#[allow(unused_imports)]
pub use numerical::optimize::*;

#[allow(unused_imports)]
pub use special::legendre::*;

#[allow(unused_imports)]
pub use redox::redoxable::*;
