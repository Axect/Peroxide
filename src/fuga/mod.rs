//! Choose what you want.
//! 
//! # Philosophy
//! 
//! Numerical algorithms are neglected in many codes. 
//! However, it is very important which algorithm is used for precise research and important numerical computation.
//! `fuga` is the best for you.
//! 
//! # Usage
//! 
//! ```ignore
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//! 
//! // Then you can use everyting in peroxide.
//! ```
//! 
//! # Compare with `prelude`
//! 
//! * Norm
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!    
//! fn main() {
//!     let a = c!(1, 2, 3);
//!     let l1 = a.norm(Norm::L1);
//!     let l2 = a.norm(Norm::L2);
//!     let l_inf = a.norm(Norm::LInf);
//!     
//!     assert_eq!(l1, 6f64);
//!     assert_eq!(l2, 14f64.sqrt());
//!     assert_eq!(l_inf, 3f64);
//! }
//! ```
//! 
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::prelude::*;
//! 
//! fn main() {
//!     let a = c!(1, 2, 3);
//!     let l2 = a.norm();      // L2 is default vector norm
//!     // prelude can't compute l1 norm, l_inf norm
//!     assert_eq!(l2, 14f64.sqrt());
//! }
//! ```
//! 
//! * Numerical integration
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

#[allow(unused_imports)]
pub use crate::macros::{
    julia_macro::*, 
    matlab_macro::*, 
    r_macro::*
};

pub use crate::traits::{
    math::{Vector, Normed, InnerProduct, Norm, LinearOp},
    fp::{FPVector, FPMatrix},
    mutable::{MutFP, MutMatrix},
    general::Algorithm,
    num::{PowOps, ExpLogOps, TrigOps, Real, Number, NumberVector},
    pointer::{Redox, Oxide, MatrixPtr},
    sugar::{VecOps, Scalable, ScalableMut},
};

#[allow(unused_imports)]
pub use crate::structure::{
    matrix::*,
    vector::*,
    dual::*,
    polynomial::*,
    hyper_dual::*,
};

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
    eigen::*,
    integral::*,
};

#[allow(unused_imports)]
pub use crate::ml::reg::*;

#[allow(unused_imports)]
#[cfg(feature = "dataframe")]
pub use crate::structure::dataframe::*;

#[allow(unused_imports)]
#[cfg(feature = "plot")]
pub use crate::util::plot::*;

// =============================================================================
// Enums
// =============================================================================
pub use crate::numerical::integral::Integral::{GaussLegendre, NewtonCotes};
pub use crate::traits::num::Number::{D, F};
pub use crate::statistics::stat::QType::{Type1, Type2, Type3, Type4, Type5, Type6, Type7, Type8, Type9};
pub use crate::structure::matrix::{
    Form::{Diagonal, Identity},
    SolveKind::{LU, WAZ},
};
