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
//! // Then you can use everything in peroxide.
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

pub use peroxide_num::{ExpLogOps, PowOps, TrigOps};

pub use crate::traits::{
    fp::{FPMatrix, FPVector},
    general::Algorithm,
    math::{InnerProduct, LinearOp, MatrixProduct, Norm, Normed, Vector, VectorProduct},
    mutable::{MutFP, MutMatrix},
    num::Real,
    pointer::{MatrixPtr, Oxide, Redox, RedoxCommon},
    stable::StableFn,
    sugar::{Scalable, ScalableMut, VecOps, ConvToMat},
};

#[allow(unused_imports)]
pub use crate::structure::{
    matrix::*, 
    polynomial::*, 
    vector::*, 
    dataframe::*,
    ad::*,
    //complex::C64,
};

pub use crate::util::{api::*, low_level::*, non_macro::*, print::*, useful::*, wrapper::*};

#[allow(unused_imports)]
pub use crate::statistics::{dist::*, ops::*, rand::*, stat::*};

#[allow(unused_imports)]
pub use crate::special::function::*;

#[allow(unused_imports)]
pub use crate::numerical::{
    eigen::*, integral::*, interp::*, ode::*, optimize::*, root::*, spline::*, utils::*,
};

#[allow(unused_imports)]
pub use crate::ml::reg::*;

#[allow(unused_imports)]
#[cfg(feature = "plot")]
pub use crate::util::plot::*;

// =============================================================================
// Enums
// =============================================================================
pub use crate::numerical::integral::Integral::{
    GaussLegendre, 
    NewtonCotes,
    G7K15,
    G10K21,
    G15K31,
    G20K41,
    G25K51,
    G30K61,
    G7K15R,
    G10K21R,
    G15K31R,
    G20K41R,
    G25K51R,
    G30K61R,
};
pub use crate::statistics::stat::QType::{
    Type1, Type2, Type3, Type4, Type5, Type6, Type7, Type8, Type9,
};
pub use crate::structure::matrix::{
    Form::{Diagonal, Identity},
    SolveKind::{LU, WAZ},
    UPLO::{Upper, Lower}
};
pub use crate::structure::dataframe::DType::*;
pub use crate::structure::ad::AD::*;
pub use crate::numerical::spline::SlopeMethod::{Akima, Quadratic};
pub use crate::statistics::stat::Metric::*;

#[cfg(feature="parquet")]
pub use arrow2::io::parquet::write::CompressionOptions;
