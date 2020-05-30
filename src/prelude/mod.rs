#[allow(unused_imports)]
pub use crate::macros::{julia_macro::*, matlab_macro::*, r_macro::*};

pub mod simpler;

pub use crate::traits::{
    math::{Vector, InnerProduct, LinearOp},
    fp::{FPVector, FPMatrix},
    mutable::{MutFP, MutMatrix},
    general::Algorithm,
    num::{PowOps, ExpLogOps, TrigOps, Real, Number, NumberVector},
    raw::RawMatrix,
};

pub use simpler::SimpleNorm;

#[allow(unused_imports)]
pub use crate::structure::{
    matrix::*,
    vector::*,
    dual::*,
    polynomial::*,
    hyper_dual::*,
};

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
pub use crate::special::{
    function::*,
};

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