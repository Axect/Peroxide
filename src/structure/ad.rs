use peroxide_ad::{
    ad_struct_def,
    ad_display,
    ad_impl,
    ad_impl_neg,
    ad_impl_add,
    ad_impl_sub,
};
use std::ops::{Neg, Add, Sub, Mul, Div};

ad_struct_def!();
ad_display!();
ad_impl!();
ad_impl_neg!();
ad_impl_add!();
ad_impl_sub!();
