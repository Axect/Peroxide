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
pub use crate::macros::{
    julia_macro::*, 
    matlab_macro::*, 
    r_macro::*
};

pub mod simpler;

pub use crate::traits::{
    math::{Vector, InnerProduct, LinearOp},
    fp::{FPVector, FPMatrix},
    mutable::{MutFP, MutMatrix},
    general::Algorithm,
    num::{PowOps, ExpLogOps, TrigOps, Real, Number, NumberVector},
    pointer::{Redox, Oxide, MatrixPtr},
};

pub use simpler::SimpleNorm;

#[allow(unused_imports)]
pub use crate::structure::{
    matrix::{
        Matrix,
        PQLU,
        WAZD,
        QR,
        Row,
        Col,
        gemv,
        gemm,
        diag,
        inv_l,
        inv_u,
        Shape,
        ml_matrix,
        py_matrix,
        r_matrix,
        combine,
        gen_householder
    },
    vector::*,
    dual::*,
    polynomial::*,
    hyper_dual::*,
};

pub use simpler::{SimplerLinearAlgebra, solve};

#[allow(unused_imports)]
pub use crate::util::{
    print::*,
    api::*,
    low_level::*,
    non_macro::*,
    useful::*,
    wrapper::*,
};

#[allow(unused_imports)]
pub use crate::statistics::{
    rand::*,
    dist::*,
    ops::*,
    stat::*,
};

#[allow(unused_imports)]
pub use crate::special::function::*;

#[allow(unused_imports)]
pub use crate::numerical::{
    ode::*,
    newton::*,
    optimize::*,
    spline::*,
    utils::*,
    interp::*,
    eigen::Eigen,
};

pub use simpler::{
    integrate,
    eigen,
};
