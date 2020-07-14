//! `peroxide` is comprehensive numerical library for Rust.
//!
//! ## Components
//!
//! `peroxide` has various components for scientific computation.
//!
//! * Linear Algebra (with BLAS & LAPACK)
//!     * [Matrix](structure/matrix/index.html) operations
//!         * `+,-,*,/`
//!         * LU, Determinant, Inverse
//!         * QR Decomposition
//!         * Reduced Row Echelon Form
//!     * [Vector](structure/vector/index.html) operations
//!     * [Eigenvalue, Eigenvector](numerical/eigen/index.html) algorithms 
//! * Statistics
//!     * [Statistical operations](statistics/stat/index.html)
//!         * `mean, var, sd`
//!         * `factorial, P, C, H`
//!     * [Distributions](statistics/dist/index.html)
//!         * Bernoulli
//!         * Uniform
//!         * Normal
//!         * Gamma
//!         * Beta
//!         * Student's t
//! * [Special functions](special/function/index.html) (Using `puruspe` crate)
//!     * Gaussian
//!     * Gamma
//!     * Beta
//!     * Error
//!     * Incomplete Gamma
//!     * Incomplete Beta
//! * Automatic Differentiation
//!     * <del>[Dual number system](structure/dual/index.html)</del>
//!     * <del>[Hyper dual number system](structure/hyper_dual/index.html)</del>
//!     * [Taylor mode forward AD](structure/ad/index.html)
//! * Numerical Utils
//!     * [Interpolation](numerical/interp/index.html)
//!     * [Spline](numerical/spline/index.html)
//!     * [Polynomial](structure/polynomial/index.html)
//!     * [Lanczos Approximation](special/lanczos/index.html)
//!     * [Numerical Integrations](numerical/integral/index.html)
//! * [Optimization](numerical/optimize/index.html)
//!     * Gradient Descent
//!     * Levenberg-Marquardt
//! * [Differential Equations](numerical/ode/index.html)
//!     * Explicit
//!         * Runge-Kutta 4th order
//!         * Euler methods
//!     * Implicit
//!         * Backward Euler
//!         * Gauss-Legendre 4th order
//! * Communication with Python
//!     * [Plot with `matplotlib`](util/plot/index.html)
//! * [DataFrame](structure/dataframe/index.html)
//!     * Read & Write with `netcdf` or `csv` format
//! * Macros
//!     * [R macros](macros/r_macro/index.html)
//!     * [Matlab macros](macros/matlab_macro/index.html)
//!     * [Julia macros](macros/julia_macro/index.html)
//!
//! And all these things are built on mathematical traits.
//! 
//! * Traits
//!     * [Functional Programming tools](traits/fp/index.html)
//!     * [General algorithms](traits/general/index.html)
//!     * [Mathematics](traits/math/index.html)
//!     * [Mutable tools](traits/mutable/index.html)
//!     * [Number & Real](traits/num/index.html)
//!     * [Pointer](traits/pointer/index.html)
//! 
//! ## Quick Start
//!
//! ### Cargo.toml
//!
//! * To use `peroxide`, you should edit `Cargo.toml`
//! * Current document version is corresponding to `0.25.0`
//!
//! 1. Default
//!     ```toml
//!     [dependencies]
//!     peroxide = "0.25"
//!     ```
//! 2. OpenBLAS & SIMD
//!     ```toml
//!     [dependencies.peroxide]
//!     version = "0.25"
//!     default-features = false
//!     features = ["O3"]
//!     ```
//! 3. Plot
//!     ```toml
//!     [dependencies.peroxide]
//!     version = "0.25"
//!     default-features = false
//!     features = ["plot"]
//!     ```
//! 4. DataFrame
//!     ```toml
//!     [dependencies.peroxide]
//!     version = "0.25"
//!     default-features = false
//!     features = ["dataframe"]
//!     ```
//! 5. Together
//!     ```toml
//!     [dependencies.peroxide]
//!     version = "0.25"
//!     default-features = false
//!     features = ["O3", "plot", "dataframe"]
//!     ```
//!
//! ## Import all at once
//!
//! Peroxide has two options.
//! 
//! * [`prelude`](prelude/inde.html) : To simple use
//! * [`fuga`](fuga/index.html) : To control numerical algorithms
//! 
//! To see differences, follow above two links.
//! 
//! You can import all functions & structures at once
//!
//! * `prelude`
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::prelude::*;
//! 
//! fn main() {
//!     // Write what you want
//! }
//! ```
//! 
//! * `fuga`
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//! 
//! fn main() {
//!     // Write what you want
//! }
//! ```
//!
//! ## Useful tips for features
//!
//! * After `0.23.0`, there are two options - `fuga`, `prelude`. Choose proper option for your computations.
//! * After `0.21.4`, if size of matrix is smaller than `1000 x 1000`, default is more effective than `O3` feature.
//! * To plot, use `dataframe` to export data as netcdf format and use python to draw plot.
//!     * `plot` feature has limited plot abilities.
//!     * There is a template of python code. - [Socialst](https://github.com/Axect/Socialst/blob/master/Templates/PyPlot_Template/nc_plot.py)

#[cfg(feature = "O3")]
extern crate blas;

#[cfg(feature = "O3")]
extern crate lapack;

#[cfg(feature = "plot")]
extern crate pyo3;

#[cfg(feature = "simd")]
extern crate packed_simd;

#[cfg(feature = "serde")]
extern crate serde;

extern crate rand;

#[cfg(feature = "dataframe")]
extern crate indexmap;

#[cfg(feature = "dataframe")]
extern crate netcdf;

#[cfg(feature = "dataframe")]
extern crate json;

extern crate order_stat;

extern crate puruspe;

extern crate matrixmultiply;

#[macro_use]
extern crate peroxide_ad;

#[macro_use]
pub mod macros;

pub mod statistics;
pub mod structure;
pub mod ml;
pub mod numerical;
pub mod special;
pub mod util;
pub mod traits;
pub mod fuga;
pub mod prelude;
