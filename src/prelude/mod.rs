//! Do not disturbed. Just use.
//!
//! # Philosophy
//!
//! For complicated computations like as physics, statistics and etc, too many options of library disturbes theory.
//! Many computations where numerical algorithms are not very critical do not require many options.
//! L2 norm is enough, and what integration algorithms you use is not important.
//! `prelude` makes you free.
//!
//! # Usage
//!
//! ```ignore
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::prelude::*;
//!
//! // Then you can use almost everything in peroxide.
//! ```
//!
//! # Compare with `fuga`
//!
//! * Norm
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::prelude::*;
//!
//! fn main() {
//!     let a = c!(1, 2, 3);
//!     let l2 = a.norm();      // L2 is default vector norm
//!
//!     assert_eq!(l2, 14f64.sqrt());
//! }
//! ```
//!
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = c!(1, 2, 3);
//!     let l2 = a.norm(Norm::L2);
//!     assert_eq!(l2, 14f64.sqrt());
//! }
//! ```
//!
//! * Numerical integration
//!
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::prelude::*;
//! use std::f64::consts::PI;
//!
//! fn main() {
//!     let sin = |x: f64| x.sin();
//!     integrate(sin, (0f64, PI)).print();
//!     // Default integration = G7K15R(1e-4, 20)
//! }
//! ```
//!
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//! use std::f64::consts::PI;
//!
//! fn main() {
//!     let sin = |x: f64| x.sin();
//!     integrate(sin, (0f64, PI), G7K15R(1e-4, 20)).print();
//! }
//! ```
//!
//! * Solve
//!
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = ml_matrix("1 2;3 4");
//!     let b = c!(3, 7);
//!     a.solve(&b, LU).print();    // [1, 1]
//!     a.solve(&b, WAZ).print();   // [1, 1]
//! }
//! ```
//!
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::prelude::*;
//!
//! fn main() {
//!     let a = ml_matrix("1 2;3 4");
//!     let b = c!(3, 7);
//!     // Prelude can only solve with LU
//!     a.solve(&b).print();    // [1, 1]
//! }
//! ```
//!
//! * DataFrame with Parquet
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let x = seq(0, 1, 0.1);
//!     let y = x.fmap(|t| t.powi(2));
//!
//!     let mut df = DataFrame::new(vec![]);
//!     df.push("x", Series::new(x));
//!     df.push("y", Series::new(y));
//!
//!     df.print();
//!
//!     # #[cfg(feature="parquet")] {
//!     df.write_parquet("example_data/test.parquet", SNAPPY).unwrap();
//!     # }
//! }
//! ```
//!
//! ```
//! extern crate peroxide;
//! use peroxide::prelude::*;
//!
//! fn main() {
//!     let x = seq(0, 1, 0.1);
//!     let y = x.fmap(|t| t.powi(2));
//!
//!     let mut df = DataFrame::new(vec![]);
//!     df.push("x", Series::new(x));
//!     df.push("y", Series::new(y));
//!
//!     df.print();
//!
//!     # #[cfg(feature="parquet")] {
//!     df.write_parquet("example_data/test.parquet").unwrap();
//!     # }
//! }
//! ```

#[allow(unused_imports)]
pub use crate::macros::{julia_macro::*, matlab_macro::*, r_macro::*};

pub use peroxide_ad::{ad_closure, ad_function};

pub mod simpler;

pub use crate::traits::{
    fp::{FPMatrix, FPVector},
    general::Algorithm,
    math::{InnerProduct, LinearOp, MatrixProduct, Vector, VectorProduct},
    matrix::{MatrixTrait, PQLU, QR, WAZD},
    mutable::{MutFP, MutMatrix},
    num::Real,
    pointer::{MatrixPtr, Oxide, Redox, RedoxCommon},
    sugar::{ConvToMat, Scalable, ScalableMut, VecOps},
};

pub use peroxide_num::{ExpLogOps, PowOps, TrigOps};

pub use simpler::SimpleNorm;

#[cfg(feature = "csv")]
pub use crate::structure::dataframe::WithCSV;
#[allow(unused_imports)]
pub use crate::structure::{
    ad::AD::*,
    ad::*,
    dataframe::{
        DType, DTypeArray, DTypeValue, DataFrame, Scalar, Series, TypedScalar, TypedVector,
    },
    matrix::{
        combine, diag, gemm, gemv, gen_householder, inv_l, inv_u, matrix, ml_matrix, py_matrix,
        r_matrix, Col, Matrix, Row, Shape,
    },
    polynomial::{lagrange_polynomial, legendre_polynomial, poly, Calculus, Polynomial},
    vector::*,
};

#[cfg(feature = "nc")]
pub use crate::structure::dataframe::WithNetCDF;

#[cfg(feature = "complex")]
#[allow(ambiguous_glob_reexports)]
#[allow(unused_imports)]
pub use crate::complex::{integral::*, matrix::*, vector::*, C64};

pub use simpler::{solve, SimplerLinearAlgebra};

#[allow(unused_imports)]
pub use crate::util::{api::*, low_level::*, non_macro::*, print::*, useful::*, wrapper::*};

#[allow(unused_imports)]
pub use crate::statistics::{dist::*, ops::*, rand::*, stat::*};

#[allow(unused_imports)]
pub use crate::special::function::{
    beta, erf, erfc, gamma, gaussian, inc_beta, inc_gamma, inv_erf, inv_erfc, inv_inc_beta,
    inv_inc_gamma, ln_gamma, phi, poch,
};

#[allow(unused_imports)]
pub use crate::numerical::{
    eigen::Eigen,
    interp::*,
    ode::*,
    optimize::*,
    root::*,
    spline::{cubic_spline, CubicHermiteSpline, CubicSpline, Spline},
    utils::*,
};

pub use simpler::{
    chebyshev_polynomial, cubic_hermite_spline, eigen, integrate, lambert_w0, lambert_wm1,
};

#[allow(unused_imports)]
pub use crate::statistics::stat::Metric::*;

#[cfg(feature = "parquet")]
pub use simpler::SimpleParquet;

#[cfg(feature = "plot")]
pub use crate::util::plot::*;

pub use anyhow;
pub use paste;
pub use rand::prelude::*;
