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
    raw::RawMatrix,
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
    eigen::*,
    integral::*,
};

#[allow(unused_imports)]
pub use crate::ml::{
    reg::*,
};

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