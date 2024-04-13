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
//!     df.write_parquet("example_data/test.parquet", CompressionOptions::Uncompressed).unwrap();
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

pub use peroxide_ad::{ad_function, ad_closure};

pub mod simpler;

pub use crate::traits::{
    fp::{FPMatrix, FPVector},
    general::Algorithm,
    math::{InnerProduct, LinearOp, MatrixProduct, Vector, VectorProduct},
    mutable::{MutFP, MutMatrix},
    num::Real,
    pointer::{MatrixPtr, Oxide, Redox, RedoxCommon},
    sugar::{Scalable, ScalableMut, VecOps, ConvToMat},
};

pub use peroxide_num::{ExpLogOps, TrigOps, PowOps};

pub use simpler::SimpleNorm;

#[allow(unused_imports)]
pub use crate::structure::{
    ad::*,
    ad::AD::*,
    matrix::{
        combine, diag, gemm, gemv, gen_householder, inv_l, inv_u, matrix, ml_matrix, py_matrix,
        r_matrix, Col, Matrix, Row, Shape, PQLU, QR, WAZD,
    },
    polynomial::{Polynomial,poly,Calculus,lagrange_polynomial,legendre_polynomial},
    vector::*,
    dataframe::{
        DataFrame, DType, DTypeArray, DTypeValue, Series, Scalar, TypedScalar, TypedVector
    },
    //complex::C64,
};
#[cfg(feature="csv")]
pub use crate::structure::dataframe::WithCSV;

#[cfg(feature="nc")]
pub use crate::structure::dataframe::WithNetCDF;

pub use simpler::{solve, SimplerLinearAlgebra};

#[allow(unused_imports)]
pub use crate::util::{api::*, low_level::*, non_macro::*, print::*, useful::*, wrapper::*};

#[allow(unused_imports)]
pub use crate::statistics::{dist::*, ops::*, rand::*, stat::*};

#[allow(unused_imports)]
pub use crate::special::function::*;

#[allow(unused_imports)]
pub use crate::numerical::{
    eigen::Eigen,
    interp::*,
    ode::*,
    optimize::*,
    root::*,
    spline::{cubic_spline, CubicSpline, CubicHermiteSpline, Spline},
    utils::*,
};

pub use simpler::{eigen, integrate, chebyshev_polynomial, cubic_hermite_spline};

#[allow(unused_imports)]
pub use crate::statistics::stat::Metric::*;

#[cfg(feature="parquet")]
pub use simpler::SimpleParquet;

#[cfg(feature="plot")]
pub use crate::util::plot::*;
