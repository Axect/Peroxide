//! `peroxide` is a comprehensive numerical computing library for Rust.
//!
//! Scientific computing workflows rarely stay inside a single numerical domain: a typical pipeline may solve an ODE, interpolate the result onto an irregular grid, evaluate a quadrature, find a root of a derived quantity, and serialize the output for downstream analysis.
//! In Python, the SciPy + NumPy stack consolidates these capabilities under one namespace.
//! In Rust the equivalent functionality is spread across specialist crates whose type systems and trait hierarchies do not always align.
//!
//! Peroxide fills that integration gap.
//! It bundles linear algebra (with optional BLAS / LAPACK acceleration), ODE integrators, quadrature, splines, root finding, optimization, statistics and distributions, const-generic forward-mode automatic differentiation (`Jet<N>`), a `DataFrame` with parquet / NetCDF / CSV I/O, and an R / MATLAB / Python style macro surface.
//! All of these are designed to interoperate through a shared `Real` trait and a compile-time-constant Butcher-tableau interface.
//!
//! The crate is aimed at researchers and engineers who want a batteries-included numerical toolbox in Rust without composing several specialist crates and reconciling their conventions.
//!
//! See the [`examples/`](https://github.com/Axect/Peroxide/tree/master/examples) directory for 40+ self-contained programs spanning every component, and the companion [Peroxide_Gallery](https://github.com/Axect/Peroxide_Gallery) repository for longer worked examples (plotting, splines, ODE applications, ...).
//!
//! ## Components
//!
//! `peroxide` has various components for scientific computation.
//!
//! - Linear Algebra (with BLAS & LAPACK)
//!   - [Matrix](structure/matrix/index.html) operations
//!     - `+,-,*,/`
//!     - LU Decomposition, Determinant, Inverse
//!     - QR Decomposition (`O3` feature needed)
//!     - Singular Value Decomposition (`O3` feature needed)
//!     - Cholesky Decomposition (`O3` feature needed)
//!     - Reduced Row Echelon Form
//!   - [Vector](structure/vector/index.html) operations
//!   - [Eigenvalue, Eigenvector](numerical/eigen/index.html) algorithms
//! - Statistics
//!   - [Statistical operations](statistics/stat/index.html)
//!     - `mean, var, sd`
//!     - `factorial, P, C, H`
//!     - Confusion Matrix & metrics
//!   - [Distributions](statistics/dist/index.html)
//!     - Bernoulli
//!     - Uniform
//!     - Binomial
//!     - Normal
//!     - Gamma
//!     - Beta
//!     - Student's t
//!     - LogNormal
//! - [Special functions](special/function/index.html) (Using `puruspe` crate)
//!   - Gaussian
//!   - Gamma
//!   - Beta
//!   - Error
//!   - Incomplete Gamma
//!   - Incomplete Beta
//!   - Lambert W
//! - Automatic Differentiation
//!   - [Const-generic `Jet<N>` forward AD](structure/ad/index.html) (arbitrary-order Taylor mode)
//!   - Type aliases: `Dual = Jet<1>`, `HyperDual = Jet<2>`
//!   - `#[ad_function]` proc macro for automatic gradient/hessian
//! - Numerical Utils
//!   - [Interpolation](numerical/interp/index.html)
//!   - [Spline](numerical/spline/index.html)
//!   - [Polynomial](structure/polynomial/index.html)
//!   - [Lanczos Approximation](special/lanczos/index.html)
//!   - [Numerical Integrations](numerical/integral/index.html)
//! - [Optimization](numerical/optimize/index.html)
//!   - Gradient Descent
//!   - Levenberg-Marquardt
//! - [Root Finding](numerical/root/index.html)
//!   - Bisection
//!   - False Position
//!   - Secant
//!   - Newton
//!   - Broyden
//! - [Ordinary Differential Equations](numerical/ode/index.html)
//!   - Explicit
//!     - Ralston's 3rd order
//!     - Runge-Kutta 4th order
//!     - Ralston's 4th order
//!     - Runge-Kutta 5th order
//!   - Embedded
//!     - Bogacki-Shampine 3(2)
//!     - Runge-Kutta-Fehlberg 4(5)
//!     - Dormand-Prince 5(4)
//!     - Tsitouras 5(4)
//!     - Runge-Kutta-Fehlberg 7(8)
//!   - Implicit
//!     - Gauss-Legendre 4th order
//! - Communication with Python
//!   - [Plot with `matplotlib`](util/plot/index.html)
//! - [DataFrame](structure/dataframe/index.html)
//!   - Read & Write with `parquet` or `netcdf` or `csv` format
//!   - Shape & info: `nrow`, `ncol`, `shape`, `dtypes`, `is_empty`, `contains`
//!   - Row operations: `head`, `tail`, `slice`
//!   - Column operations: `select`, `rename`, `column_names`, `select_dtypes`
//!   - Statistics: `describe`, `sum`, `mean`
//!   - Series statistics: `sum`, `mean`, `var`, `sd`, `min`, `max`
//! - Macros
//!   - [R macros](macros/r_macro/index.html)
//!   - [Matlab macros](macros/matlab_macro/index.html)
//!   - [Julia macros](macros/julia_macro/index.html)
//!
//! And all these things are built on mathematical traits.
//!
//! - Traits
//!   - [Functional Programming tools](traits/fp/index.html)
//!   - [General algorithms](traits/general/index.html)
//!   - [Mathematics](traits/math/index.html)
//!   - [Mutable tools](traits/mutable/index.html)
//!   - [Number & Real](traits/num/index.html)
//!   - [Pointer](traits/pointer/index.html)
//!   - [Stable](traits/stable/index.html)
//! ## Import all at once
//!
//! Peroxide has two options.
//!
//! * [`prelude`](prelude/index.html) : To simple use
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
//! * If you want to use _QR_, _SVD_, or _Cholesky Decomposition_, you should use the `O3` feature. These decompositions are not implemented in the `default` feature.
//!
//! * If you want to save your numerical results, consider using the `parquet` or `nc` features, which correspond to the `parquet` and `netcdf` file formats, respectively. These formats are much more efficient than `csv` and `json`.
//!
//! * For plotting, it is recommended to use the `plot` feature. However, if you require more customization, you can use the `parquet` or `nc` feature to export your data in the parquet or netcdf format and then use Python to create the plots.
//!
//!     * To read parquet files in Python, you can use the `pandas` and `pyarrow` libraries.
//!
//!     * A template for Python code that works with netcdf files can be found in the [Socialst](https://github.com/Axect/Socialst/blob/master/Templates/PyPlot_Template/nc_plot.py) repository.

//!

#[cfg(feature = "O3")]
extern crate blas;

#[cfg(feature = "O3")]
extern crate lapack;

// Force link-directive propagation when a backend convenience feature
// (`O3-openblas` / `O3-accelerate` / `O3-mkl` / `O3-netlib`) activates
// `blas-src` / `lapack-src` via Cargo's optional-dep feature shim.
#[cfg(feature = "blas-src")]
extern crate blas_src as _;
#[cfg(feature = "lapack-src")]
extern crate lapack_src as _;

#[cfg(feature = "plot")]
extern crate pyo3;

#[cfg(feature = "serde")]
extern crate serde;

extern crate rand;

// extern crate json;

extern crate order_stat;

extern crate puruspe;

extern crate matrixmultiply;

#[cfg(feature = "nc")]
extern crate netcdf;

extern crate peroxide_ad;

#[macro_use]
pub mod macros;

pub mod fuga;
pub mod ml;
pub mod numerical;
pub mod prelude;
pub mod special;
pub mod statistics;
pub mod structure;
pub mod traits;
pub mod util;

#[cfg(feature = "complex")]
pub mod complex;

#[cfg(feature = "parallel")]
extern crate rayon;
