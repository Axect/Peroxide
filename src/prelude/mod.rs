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
//!     // Default integration = GaussLegendre(15)
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
//!     integrate(sin, (0f64, PI), GaussLegendre(15)).print();
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

#[allow(unused_imports)]
pub use crate::macros::{julia_macro::*, matlab_macro::*, r_macro::*};

pub mod simpler;

pub use crate::traits::{
    fp::{FPMatrix, FPVector},
    general::Algorithm,
    math::{InnerProduct, LinearOp, MatrixProduct, Vector, VectorProduct},
    mutable::{MutFP, MutMatrix},
    num::{ExpLogOps, PowOps, Real, TrigOps},
    pointer::{MatrixPtr, Oxide, Redox, RedoxCommon},
    sugar::{Scalable, ScalableMut, VecOps, ConvToMat},
};

pub use simpler::SimpleNorm;

#[allow(unused_imports)]
pub use crate::structure::{
    ad::*,
    ad::AD::*,
    matrix::{
        combine, diag, gemm, gemv, gen_householder, inv_l, inv_u, matrix, ml_matrix, py_matrix,
        r_matrix, Col, Matrix, Row, Shape, PQLU, QR, WAZD,
    },
    polynomial::*,
    vector::*,
    dataframe::*,
    complex::C64,
};

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
    root::{bisection, false_position, newton, secant},
    spline::*,
    utils::*,
};

pub use simpler::{eigen, integrate};
