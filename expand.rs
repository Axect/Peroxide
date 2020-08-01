pub mod ad {
    //! Taylor mode forward automatic differentiation
    //!
    //! **Caution**: In this documentation, `{i}` means number `i` (ex: AD{2} means AD2)
    //!
    //! ## Important Features
    //!
    //! * Can automatic differentiate up to 10th order (`AD1` ~ `AD10`)
    //! * All `AD{i}` are in stack (Guarantee high performance)
    //! * You can see `AD{i}` via `.print()`
    //! * You can use n-th derivative value of `AD{i}` via `.d{i}`
    //!
    //! ## Implemented Traits for `AD{i}`
    //!
    //! * `#[derive(Debug, Copy, Clone, PartialEq, Default)]`
    //! * `std::fmt::Display`
    //! * `From<AD{j}>`, `From<&'a AD{j}>`
    //! * `IntoIterator<Item = f64>`
    //! * `FromIterator<f64>`
    //! * `Index`, `IndexMut`
    //! * `std::ops::{Neg, Add, Sub, Mul, Div}`
    //! * `peroxide::traits::num::{PowOps, ExpLogOps, TrigOps}`
    //!
    //! ## Iterator of `AD{i}`
    //!
    //! There are three iterators.
    //!
    //! * `ADIntoIter{i}`
    //! * `ADIter{i}<'a>`
    //! * `ADIterMut{i}<'a>`
    //!
    //! Each implements `DoubleEndedIterator`, `ExactSizeIterator` also.
    //!
    //! ## Methods
    //!
    //! * `new(d0: f64, ... , d{i}: f64) -> AD{i}`
    //! * `print(&self)`
    //! * `iter(&self) -> ADIter{i}`
    //! * `iter_mut(&self) -> ADIterMut{i}`
    //! * `len(&self) -> usize`
    //!
    //! ## Implemented Operations
    //!
    //! * `Add, Sub, Mul, Div`
    //! * `sin, cos, tan`
    //! * `sinh, cosh, tanh`
    //! * `sin_cos`, `sinh_cosh`
    //! * `exp, ln, log, log2, log10`
    //! * `powi, powf, sqrt`
    //!
    //! ## Not yet implemented
    //!
    //! * `asin`, `acos`, `atan`
    //! * `asinh`, `acosh`, `atanh`
    //! * `pow`
    //!
    //! ## Usage
    //!
    //! ### Construction
    //!
    //! ```
    //! extern crate peroxide;
    //! use peroxide::fuga::*;
    //!
    //! fn main() {
    //!     // Declare x where x = 2
    //!     let a = AD1::new(2f64, 1f64);
    //!     // Declare x^2 where x = 2
    //!     let b = AD2::new(4f64, 4f64, 2f64);
    //!     // Convert AD1 -> AD2
    //!     let c = AD2::from(a);
    //!     // Zeros
    //!     let d = AD2::default();
    //!
    //!     assert_eq!(c, AD2::new(2f64, 1f64, 0f64));
    //!     assert_eq!(d, AD2::new(0f64, 0f64, 0f64));
    //! }
    //! ```
    //!
    //! ### Operation
    //!
    //! For every binary operation, it returns higher order AD
    //! (ex: AD1 + AD2 = AD2)
    //!
    //! ```
    //! extern crate peroxide;
    //! use peroxide::fuga::*;
    //!
    //! fn main() {
    //!     let a = AD1::new(2f64, 1f64);       // x        at x = 2
    //!     let b = AD2::new(4f64, 4f64, 2f64); // x^2      at x = 2
    //!     let c = a + b;                      // x^2 + x  at x = 2
    //!     let d = a * b;                      // x^3      at x = 2
    //!     let e = a / b;                      // 1/x      at x = 2
    //!     assert_eq!(c, AD2::new(6f64, 5f64, 2f64));
    //!     assert_eq!(d, AD2::new(8f64, 12f64, 12f64));
    //!     assert_eq!(e, AD2::new(0.5, -0.25, 0.25));
    //! }
    //! ```
    //!
    //! ### Generic
    //!
    //! * All of `AD{i}` implements `AD` trait
    //!
    //! ```
    //! extern crate peroxide;
    //! use peroxide::fuga::*;
    //!
    //! fn main() {
    //!     let a = AD1::new(2f64, 1f64);
    //!     let b = AD2::new(4f64, 4f64, 2f64);
    //!     assert_eq!(f(a, b), AD1::new(6f64, 5f64));
    //! }
    //!
    //! fn f<T: AD, S: AD>(x: T, y: S) -> T {
    //!     T::from(x.to_ad2() + y.to_ad2())
    //! }
    //! ```
    use crate::statistics::ops::C;
    use crate::traits::{
        num::{ExpLogOps, PowOps, TrigOps},
        stable::StableFn,
    };
    use peroxide_ad::{
        ad_display, ad_impl, ad_impl_ad, ad_impl_add, ad_impl_div, ad_impl_double_ended_iter,
        ad_impl_exact_size_iter, ad_impl_explogops, ad_impl_from, ad_impl_from_iter, ad_impl_index,
        ad_impl_into_iter, ad_impl_iter, ad_impl_mul, ad_impl_neg, ad_impl_powops, ad_impl_sub,
        ad_impl_trigops, ad_iter_def, ad_struct_def, ad_impl_from_type, ad_impl_add_f64,
        ad_impl_sub_f64, ad_impl_mul_f64, ad_impl_div_f64, f64_impl_add_ad, f64_impl_sub_ad,
        f64_impl_mul_ad, f64_impl_div_ad, f64_impl_from_ad, ad_impl_stable_fn, def_ad,
    };
    use std::iter::FromIterator;
    use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};
    use std::marker::PhantomData;
    pub struct AD1 {
        pub d0: f64,
        pub d1: f64,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for AD1 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                AD1 {
                    d0: ref __self_0_0,
                    d1: ref __self_0_1,
                } => {
                    let mut debug_trait_builder = f.debug_struct("AD1");
                    let _ = debug_trait_builder.field("d0", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("d1", &&(*__self_0_1));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::marker::Copy for AD1 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::clone::Clone for AD1 {
        #[inline]
        fn clone(&self) -> AD1 {
            {
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                *self
            }
        }
    }
    impl ::core::marker::StructuralPartialEq for AD1 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::cmp::PartialEq for AD1 {
        #[inline]
        fn eq(&self, other: &AD1) -> bool {
            match *other {
                AD1 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                } => match *self {
                    AD1 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                    } => (*__self_0_0) == (*__self_1_0) && (*__self_0_1) == (*__self_1_1),
                },
            }
        }
        #[inline]
        fn ne(&self, other: &AD1) -> bool {
            match *other {
                AD1 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                } => match *self {
                    AD1 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                    } => (*__self_0_0) != (*__self_1_0) || (*__self_0_1) != (*__self_1_1),
                },
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::default::Default for AD1 {
        #[inline]
        fn default() -> AD1 {
            AD1 {
                d0: ::core::default::Default::default(),
                d1: ::core::default::Default::default(),
            }
        }
    }
    pub struct AD2 {
        pub d0: f64,
        pub d1: f64,
        pub d2: f64,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for AD2 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                AD2 {
                    d0: ref __self_0_0,
                    d1: ref __self_0_1,
                    d2: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("AD2");
                    let _ = debug_trait_builder.field("d0", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("d1", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("d2", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::marker::Copy for AD2 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::clone::Clone for AD2 {
        #[inline]
        fn clone(&self) -> AD2 {
            {
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                *self
            }
        }
    }
    impl ::core::marker::StructuralPartialEq for AD2 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::cmp::PartialEq for AD2 {
        #[inline]
        fn eq(&self, other: &AD2) -> bool {
            match *other {
                AD2 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                    d2: ref __self_1_2,
                } => match *self {
                    AD2 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                        d2: ref __self_0_2,
                    } => {
                        (*__self_0_0) == (*__self_1_0)
                            && (*__self_0_1) == (*__self_1_1)
                            && (*__self_0_2) == (*__self_1_2)
                    }
                },
            }
        }
        #[inline]
        fn ne(&self, other: &AD2) -> bool {
            match *other {
                AD2 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                    d2: ref __self_1_2,
                } => match *self {
                    AD2 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                        d2: ref __self_0_2,
                    } => {
                        (*__self_0_0) != (*__self_1_0)
                            || (*__self_0_1) != (*__self_1_1)
                            || (*__self_0_2) != (*__self_1_2)
                    }
                },
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::default::Default for AD2 {
        #[inline]
        fn default() -> AD2 {
            AD2 {
                d0: ::core::default::Default::default(),
                d1: ::core::default::Default::default(),
                d2: ::core::default::Default::default(),
            }
        }
    }
    pub struct AD3 {
        pub d0: f64,
        pub d1: f64,
        pub d2: f64,
        pub d3: f64,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for AD3 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                AD3 {
                    d0: ref __self_0_0,
                    d1: ref __self_0_1,
                    d2: ref __self_0_2,
                    d3: ref __self_0_3,
                } => {
                    let mut debug_trait_builder = f.debug_struct("AD3");
                    let _ = debug_trait_builder.field("d0", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("d1", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("d2", &&(*__self_0_2));
                    let _ = debug_trait_builder.field("d3", &&(*__self_0_3));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::marker::Copy for AD3 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::clone::Clone for AD3 {
        #[inline]
        fn clone(&self) -> AD3 {
            {
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                *self
            }
        }
    }
    impl ::core::marker::StructuralPartialEq for AD3 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::cmp::PartialEq for AD3 {
        #[inline]
        fn eq(&self, other: &AD3) -> bool {
            match *other {
                AD3 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                    d2: ref __self_1_2,
                    d3: ref __self_1_3,
                } => match *self {
                    AD3 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                        d2: ref __self_0_2,
                        d3: ref __self_0_3,
                    } => {
                        (*__self_0_0) == (*__self_1_0)
                            && (*__self_0_1) == (*__self_1_1)
                            && (*__self_0_2) == (*__self_1_2)
                            && (*__self_0_3) == (*__self_1_3)
                    }
                },
            }
        }
        #[inline]
        fn ne(&self, other: &AD3) -> bool {
            match *other {
                AD3 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                    d2: ref __self_1_2,
                    d3: ref __self_1_3,
                } => match *self {
                    AD3 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                        d2: ref __self_0_2,
                        d3: ref __self_0_3,
                    } => {
                        (*__self_0_0) != (*__self_1_0)
                            || (*__self_0_1) != (*__self_1_1)
                            || (*__self_0_2) != (*__self_1_2)
                            || (*__self_0_3) != (*__self_1_3)
                    }
                },
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::default::Default for AD3 {
        #[inline]
        fn default() -> AD3 {
            AD3 {
                d0: ::core::default::Default::default(),
                d1: ::core::default::Default::default(),
                d2: ::core::default::Default::default(),
                d3: ::core::default::Default::default(),
            }
        }
    }
    pub struct AD4 {
        pub d0: f64,
        pub d1: f64,
        pub d2: f64,
        pub d3: f64,
        pub d4: f64,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for AD4 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                AD4 {
                    d0: ref __self_0_0,
                    d1: ref __self_0_1,
                    d2: ref __self_0_2,
                    d3: ref __self_0_3,
                    d4: ref __self_0_4,
                } => {
                    let mut debug_trait_builder = f.debug_struct("AD4");
                    let _ = debug_trait_builder.field("d0", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("d1", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("d2", &&(*__self_0_2));
                    let _ = debug_trait_builder.field("d3", &&(*__self_0_3));
                    let _ = debug_trait_builder.field("d4", &&(*__self_0_4));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::marker::Copy for AD4 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::clone::Clone for AD4 {
        #[inline]
        fn clone(&self) -> AD4 {
            {
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                *self
            }
        }
    }
    impl ::core::marker::StructuralPartialEq for AD4 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::cmp::PartialEq for AD4 {
        #[inline]
        fn eq(&self, other: &AD4) -> bool {
            match *other {
                AD4 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                    d2: ref __self_1_2,
                    d3: ref __self_1_3,
                    d4: ref __self_1_4,
                } => match *self {
                    AD4 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                        d2: ref __self_0_2,
                        d3: ref __self_0_3,
                        d4: ref __self_0_4,
                    } => {
                        (*__self_0_0) == (*__self_1_0)
                            && (*__self_0_1) == (*__self_1_1)
                            && (*__self_0_2) == (*__self_1_2)
                            && (*__self_0_3) == (*__self_1_3)
                            && (*__self_0_4) == (*__self_1_4)
                    }
                },
            }
        }
        #[inline]
        fn ne(&self, other: &AD4) -> bool {
            match *other {
                AD4 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                    d2: ref __self_1_2,
                    d3: ref __self_1_3,
                    d4: ref __self_1_4,
                } => match *self {
                    AD4 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                        d2: ref __self_0_2,
                        d3: ref __self_0_3,
                        d4: ref __self_0_4,
                    } => {
                        (*__self_0_0) != (*__self_1_0)
                            || (*__self_0_1) != (*__self_1_1)
                            || (*__self_0_2) != (*__self_1_2)
                            || (*__self_0_3) != (*__self_1_3)
                            || (*__self_0_4) != (*__self_1_4)
                    }
                },
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::default::Default for AD4 {
        #[inline]
        fn default() -> AD4 {
            AD4 {
                d0: ::core::default::Default::default(),
                d1: ::core::default::Default::default(),
                d2: ::core::default::Default::default(),
                d3: ::core::default::Default::default(),
                d4: ::core::default::Default::default(),
            }
        }
    }
    pub struct AD5 {
        pub d0: f64,
        pub d1: f64,
        pub d2: f64,
        pub d3: f64,
        pub d4: f64,
        pub d5: f64,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for AD5 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                AD5 {
                    d0: ref __self_0_0,
                    d1: ref __self_0_1,
                    d2: ref __self_0_2,
                    d3: ref __self_0_3,
                    d4: ref __self_0_4,
                    d5: ref __self_0_5,
                } => {
                    let mut debug_trait_builder = f.debug_struct("AD5");
                    let _ = debug_trait_builder.field("d0", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("d1", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("d2", &&(*__self_0_2));
                    let _ = debug_trait_builder.field("d3", &&(*__self_0_3));
                    let _ = debug_trait_builder.field("d4", &&(*__self_0_4));
                    let _ = debug_trait_builder.field("d5", &&(*__self_0_5));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::marker::Copy for AD5 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::clone::Clone for AD5 {
        #[inline]
        fn clone(&self) -> AD5 {
            {
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                let _: ::core::clone::AssertParamIsClone<f64>;
                *self
            }
        }
    }
    impl ::core::marker::StructuralPartialEq for AD5 {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::cmp::PartialEq for AD5 {
        #[inline]
        fn eq(&self, other: &AD5) -> bool {
            match *other {
                AD5 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                    d2: ref __self_1_2,
                    d3: ref __self_1_3,
                    d4: ref __self_1_4,
                    d5: ref __self_1_5,
                } => match *self {
                    AD5 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                        d2: ref __self_0_2,
                        d3: ref __self_0_3,
                        d4: ref __self_0_4,
                        d5: ref __self_0_5,
                    } => {
                        (*__self_0_0) == (*__self_1_0)
                            && (*__self_0_1) == (*__self_1_1)
                            && (*__self_0_2) == (*__self_1_2)
                            && (*__self_0_3) == (*__self_1_3)
                            && (*__self_0_4) == (*__self_1_4)
                            && (*__self_0_5) == (*__self_1_5)
                    }
                },
            }
        }
        #[inline]
        fn ne(&self, other: &AD5) -> bool {
            match *other {
                AD5 {
                    d0: ref __self_1_0,
                    d1: ref __self_1_1,
                    d2: ref __self_1_2,
                    d3: ref __self_1_3,
                    d4: ref __self_1_4,
                    d5: ref __self_1_5,
                } => match *self {
                    AD5 {
                        d0: ref __self_0_0,
                        d1: ref __self_0_1,
                        d2: ref __self_0_2,
                        d3: ref __self_0_3,
                        d4: ref __self_0_4,
                        d5: ref __self_0_5,
                    } => {
                        (*__self_0_0) != (*__self_1_0)
                            || (*__self_0_1) != (*__self_1_1)
                            || (*__self_0_2) != (*__self_1_2)
                            || (*__self_0_3) != (*__self_1_3)
                            || (*__self_0_4) != (*__self_1_4)
                            || (*__self_0_5) != (*__self_1_5)
                    }
                },
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::default::Default for AD5 {
        #[inline]
        fn default() -> AD5 {
            AD5 {
                d0: ::core::default::Default::default(),
                d1: ::core::default::Default::default(),
                d2: ::core::default::Default::default(),
                d3: ::core::default::Default::default(),
                d4: ::core::default::Default::default(),
                d5: ::core::default::Default::default(),
            }
        }
    }
    impl std::fmt::Display for AD1 {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            let mut s = {
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["AD\n  d0: ", "\n"],
                    &match (&self.d0,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            };
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d1: ", "\n"],
                    &match (&self.d1,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &[""],
                &match (&s,) {
                    (arg0,) => [::core::fmt::ArgumentV1::new(
                        arg0,
                        ::core::fmt::Display::fmt,
                    )],
                },
            ))
        }
    }
    impl std::fmt::Display for AD2 {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            let mut s = {
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["AD\n  d0: ", "\n"],
                    &match (&self.d0,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            };
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d1: ", "\n"],
                    &match (&self.d1,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d2: ", "\n"],
                    &match (&self.d2,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &[""],
                &match (&s,) {
                    (arg0,) => [::core::fmt::ArgumentV1::new(
                        arg0,
                        ::core::fmt::Display::fmt,
                    )],
                },
            ))
        }
    }
    impl std::fmt::Display for AD3 {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            let mut s = {
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["AD\n  d0: ", "\n"],
                    &match (&self.d0,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            };
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d1: ", "\n"],
                    &match (&self.d1,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d2: ", "\n"],
                    &match (&self.d2,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d3: ", "\n"],
                    &match (&self.d3,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &[""],
                &match (&s,) {
                    (arg0,) => [::core::fmt::ArgumentV1::new(
                        arg0,
                        ::core::fmt::Display::fmt,
                    )],
                },
            ))
        }
    }
    impl std::fmt::Display for AD4 {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            let mut s = {
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["AD\n  d0: ", "\n"],
                    &match (&self.d0,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            };
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d1: ", "\n"],
                    &match (&self.d1,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d2: ", "\n"],
                    &match (&self.d2,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d3: ", "\n"],
                    &match (&self.d3,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d4: ", "\n"],
                    &match (&self.d4,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &[""],
                &match (&s,) {
                    (arg0,) => [::core::fmt::ArgumentV1::new(
                        arg0,
                        ::core::fmt::Display::fmt,
                    )],
                },
            ))
        }
    }
    impl std::fmt::Display for AD5 {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            let mut s = {
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["AD\n  d0: ", "\n"],
                    &match (&self.d0,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            };
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d1: ", "\n"],
                    &match (&self.d1,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d2: ", "\n"],
                    &match (&self.d2,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d3: ", "\n"],
                    &match (&self.d3,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d4: ", "\n"],
                    &match (&self.d4,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            s.push_str(&{
                let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                    &["  d5: ", "\n"],
                    &match (&self.d5,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
                res
            });
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &[""],
                &match (&s,) {
                    (arg0,) => [::core::fmt::ArgumentV1::new(
                        arg0,
                        ::core::fmt::Display::fmt,
                    )],
                },
            ))
        }
    }
    impl AD1 {
        pub fn new(d0: f64, d1: f64) -> Self {
            Self { d0, d1 }
        }
        pub fn print(&self) {
            {
                ::std::io::_print(::core::fmt::Arguments::new_v1(
                    &["", "\n"],
                    &match (&self,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
            };
        }
        pub fn iter(&self) -> ADIter1 {
            self.into_iter()
        }
        pub fn iter_mut(&mut self) -> ADIterMut1 {
            self.into_iter()
        }
        pub fn len(&self) -> usize {
            2
        }
    }
    impl AD2 {
        pub fn new(d0: f64, d1: f64, d2: f64) -> Self {
            Self { d0, d1, d2 }
        }
        pub fn print(&self) {
            {
                ::std::io::_print(::core::fmt::Arguments::new_v1(
                    &["", "\n"],
                    &match (&self,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
            };
        }
        pub fn iter(&self) -> ADIter2 {
            self.into_iter()
        }
        pub fn iter_mut(&mut self) -> ADIterMut2 {
            self.into_iter()
        }
        pub fn len(&self) -> usize {
            3
        }
    }
    impl AD3 {
        pub fn new(d0: f64, d1: f64, d2: f64, d3: f64) -> Self {
            Self { d0, d1, d2, d3 }
        }
        pub fn print(&self) {
            {
                ::std::io::_print(::core::fmt::Arguments::new_v1(
                    &["", "\n"],
                    &match (&self,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
            };
        }
        pub fn iter(&self) -> ADIter3 {
            self.into_iter()
        }
        pub fn iter_mut(&mut self) -> ADIterMut3 {
            self.into_iter()
        }
        pub fn len(&self) -> usize {
            4
        }
    }
    impl AD4 {
        pub fn new(d0: f64, d1: f64, d2: f64, d3: f64, d4: f64) -> Self {
            Self { d0, d1, d2, d3, d4 }
        }
        pub fn print(&self) {
            {
                ::std::io::_print(::core::fmt::Arguments::new_v1(
                    &["", "\n"],
                    &match (&self,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
            };
        }
        pub fn iter(&self) -> ADIter4 {
            self.into_iter()
        }
        pub fn iter_mut(&mut self) -> ADIterMut4 {
            self.into_iter()
        }
        pub fn len(&self) -> usize {
            5
        }
    }
    impl AD5 {
        pub fn new(d0: f64, d1: f64, d2: f64, d3: f64, d4: f64, d5: f64) -> Self {
            Self {
                d0,
                d1,
                d2,
                d3,
                d4,
                d5,
            }
        }
        pub fn print(&self) {
            {
                ::std::io::_print(::core::fmt::Arguments::new_v1(
                    &["", "\n"],
                    &match (&self,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                ));
            };
        }
        pub fn iter(&self) -> ADIter5 {
            self.into_iter()
        }
        pub fn iter_mut(&mut self) -> ADIterMut5 {
            self.into_iter()
        }
        pub fn len(&self) -> usize {
            6
        }
    }
    impl<'a> From<&'a AD1> for AD1 {
        fn from(high: &'a AD1) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl From<AD1> for AD2 {
        fn from(low: AD1) -> Self {
            Self::new(low.d0, low.d1, 0f64)
        }
    }
    impl<'a> From<&'a AD1> for AD2 {
        fn from(low: &'a AD1) -> Self {
            Self::new(low.d0, low.d1, 0f64)
        }
    }
    impl From<AD1> for AD3 {
        fn from(low: AD1) -> Self {
            Self::new(low.d0, low.d1, 0f64, 0f64)
        }
    }
    impl<'a> From<&'a AD1> for AD3 {
        fn from(low: &'a AD1) -> Self {
            Self::new(low.d0, low.d1, 0f64, 0f64)
        }
    }
    impl From<AD1> for AD4 {
        fn from(low: AD1) -> Self {
            Self::new(low.d0, low.d1, 0f64, 0f64, 0f64)
        }
    }
    impl<'a> From<&'a AD1> for AD4 {
        fn from(low: &'a AD1) -> Self {
            Self::new(low.d0, low.d1, 0f64, 0f64, 0f64)
        }
    }
    impl From<AD1> for AD5 {
        fn from(low: AD1) -> Self {
            Self::new(low.d0, low.d1, 0f64, 0f64, 0f64, 0f64)
        }
    }
    impl<'a> From<&'a AD1> for AD5 {
        fn from(low: &'a AD1) -> Self {
            Self::new(low.d0, low.d1, 0f64, 0f64, 0f64, 0f64)
        }
    }
    impl From<AD2> for AD1 {
        fn from(high: AD2) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl<'a> From<&'a AD2> for AD1 {
        fn from(high: &'a AD2) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl<'a> From<&'a AD2> for AD2 {
        fn from(high: &'a AD2) -> Self {
            Self::new(high.d0, high.d1, high.d2)
        }
    }
    impl From<AD2> for AD3 {
        fn from(low: AD2) -> Self {
            Self::new(low.d0, low.d1, low.d2, 0f64)
        }
    }
    impl<'a> From<&'a AD2> for AD3 {
        fn from(low: &'a AD2) -> Self {
            Self::new(low.d0, low.d1, low.d2, 0f64)
        }
    }
    impl From<AD2> for AD4 {
        fn from(low: AD2) -> Self {
            Self::new(low.d0, low.d1, low.d2, 0f64, 0f64)
        }
    }
    impl<'a> From<&'a AD2> for AD4 {
        fn from(low: &'a AD2) -> Self {
            Self::new(low.d0, low.d1, low.d2, 0f64, 0f64)
        }
    }
    impl From<AD2> for AD5 {
        fn from(low: AD2) -> Self {
            Self::new(low.d0, low.d1, low.d2, 0f64, 0f64, 0f64)
        }
    }
    impl<'a> From<&'a AD2> for AD5 {
        fn from(low: &'a AD2) -> Self {
            Self::new(low.d0, low.d1, low.d2, 0f64, 0f64, 0f64)
        }
    }
    impl From<AD3> for AD1 {
        fn from(high: AD3) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl<'a> From<&'a AD3> for AD1 {
        fn from(high: &'a AD3) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl From<AD3> for AD2 {
        fn from(high: AD3) -> Self {
            Self::new(high.d0, high.d1, high.d2)
        }
    }
    impl<'a> From<&'a AD3> for AD2 {
        fn from(high: &'a AD3) -> Self {
            Self::new(high.d0, high.d1, high.d2)
        }
    }
    impl<'a> From<&'a AD3> for AD3 {
        fn from(high: &'a AD3) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3)
        }
    }
    impl From<AD3> for AD4 {
        fn from(low: AD3) -> Self {
            Self::new(low.d0, low.d1, low.d2, low.d3, 0f64)
        }
    }
    impl<'a> From<&'a AD3> for AD4 {
        fn from(low: &'a AD3) -> Self {
            Self::new(low.d0, low.d1, low.d2, low.d3, 0f64)
        }
    }
    impl From<AD3> for AD5 {
        fn from(low: AD3) -> Self {
            Self::new(low.d0, low.d1, low.d2, low.d3, 0f64, 0f64)
        }
    }
    impl<'a> From<&'a AD3> for AD5 {
        fn from(low: &'a AD3) -> Self {
            Self::new(low.d0, low.d1, low.d2, low.d3, 0f64, 0f64)
        }
    }
    impl From<AD4> for AD1 {
        fn from(high: AD4) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl<'a> From<&'a AD4> for AD1 {
        fn from(high: &'a AD4) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl From<AD4> for AD2 {
        fn from(high: AD4) -> Self {
            Self::new(high.d0, high.d1, high.d2)
        }
    }
    impl<'a> From<&'a AD4> for AD2 {
        fn from(high: &'a AD4) -> Self {
            Self::new(high.d0, high.d1, high.d2)
        }
    }
    impl From<AD4> for AD3 {
        fn from(high: AD4) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3)
        }
    }
    impl<'a> From<&'a AD4> for AD3 {
        fn from(high: &'a AD4) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3)
        }
    }
    impl<'a> From<&'a AD4> for AD4 {
        fn from(high: &'a AD4) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3, high.d4)
        }
    }
    impl From<AD4> for AD5 {
        fn from(low: AD4) -> Self {
            Self::new(low.d0, low.d1, low.d2, low.d3, low.d4, 0f64)
        }
    }
    impl<'a> From<&'a AD4> for AD5 {
        fn from(low: &'a AD4) -> Self {
            Self::new(low.d0, low.d1, low.d2, low.d3, low.d4, 0f64)
        }
    }
    impl From<AD5> for AD1 {
        fn from(high: AD5) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl<'a> From<&'a AD5> for AD1 {
        fn from(high: &'a AD5) -> Self {
            Self::new(high.d0, high.d1)
        }
    }
    impl From<AD5> for AD2 {
        fn from(high: AD5) -> Self {
            Self::new(high.d0, high.d1, high.d2)
        }
    }
    impl<'a> From<&'a AD5> for AD2 {
        fn from(high: &'a AD5) -> Self {
            Self::new(high.d0, high.d1, high.d2)
        }
    }
    impl From<AD5> for AD3 {
        fn from(high: AD5) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3)
        }
    }
    impl<'a> From<&'a AD5> for AD3 {
        fn from(high: &'a AD5) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3)
        }
    }
    impl From<AD5> for AD4 {
        fn from(high: AD5) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3, high.d4)
        }
    }
    impl<'a> From<&'a AD5> for AD4 {
        fn from(high: &'a AD5) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3, high.d4)
        }
    }
    impl<'a> From<&'a AD5> for AD5 {
        fn from(high: &'a AD5) -> Self {
            Self::new(high.d0, high.d1, high.d2, high.d3, high.d4, high.d5)
        }
    }
    pub struct ADIntoIter1 {
        ad: AD1,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for ADIntoIter1 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIntoIter1 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIntoIter1");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIter1<'a> {
        ad: &'a AD1,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIter1<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIter1 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIter1");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIterMut1<'a> {
        ad: &'a mut AD1,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIterMut1<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIterMut1 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIterMut1");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIntoIter2 {
        ad: AD2,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for ADIntoIter2 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIntoIter2 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIntoIter2");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIter2<'a> {
        ad: &'a AD2,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIter2<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIter2 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIter2");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIterMut2<'a> {
        ad: &'a mut AD2,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIterMut2<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIterMut2 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIterMut2");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIntoIter3 {
        ad: AD3,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for ADIntoIter3 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIntoIter3 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIntoIter3");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIter3<'a> {
        ad: &'a AD3,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIter3<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIter3 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIter3");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIterMut3<'a> {
        ad: &'a mut AD3,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIterMut3<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIterMut3 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIterMut3");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIntoIter4 {
        ad: AD4,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for ADIntoIter4 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIntoIter4 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIntoIter4");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIter4<'a> {
        ad: &'a AD4,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIter4<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIter4 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIter4");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIterMut4<'a> {
        ad: &'a mut AD4,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIterMut4<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIterMut4 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIterMut4");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIntoIter5 {
        ad: AD5,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::fmt::Debug for ADIntoIter5 {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIntoIter5 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIntoIter5");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIter5<'a> {
        ad: &'a AD5,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIter5<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIter5 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIter5");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    pub struct ADIterMut5<'a> {
        ad: &'a mut AD5,
        index: usize,
        r_index: usize,
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl<'a> ::core::fmt::Debug for ADIterMut5<'a> {
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match *self {
                ADIterMut5 {
                    ad: ref __self_0_0,
                    index: ref __self_0_1,
                    r_index: ref __self_0_2,
                } => {
                    let mut debug_trait_builder = f.debug_struct("ADIterMut5");
                    let _ = debug_trait_builder.field("ad", &&(*__self_0_0));
                    let _ = debug_trait_builder.field("index", &&(*__self_0_1));
                    let _ = debug_trait_builder.field("r_index", &&(*__self_0_2));
                    debug_trait_builder.finish()
                }
            }
        }
    }
    impl IntoIterator for AD1 {
        type Item = f64;
        type IntoIter = ADIntoIter1;
        fn into_iter(self) -> Self::IntoIter {
            ADIntoIter1 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a AD1 {
        type Item = f64;
        type IntoIter = ADIter1<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIter1 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a mut AD1 {
        type Item = f64;
        type IntoIter = ADIterMut1<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIterMut1 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl IntoIterator for AD2 {
        type Item = f64;
        type IntoIter = ADIntoIter2;
        fn into_iter(self) -> Self::IntoIter {
            ADIntoIter2 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a AD2 {
        type Item = f64;
        type IntoIter = ADIter2<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIter2 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a mut AD2 {
        type Item = f64;
        type IntoIter = ADIterMut2<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIterMut2 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl IntoIterator for AD3 {
        type Item = f64;
        type IntoIter = ADIntoIter3;
        fn into_iter(self) -> Self::IntoIter {
            ADIntoIter3 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a AD3 {
        type Item = f64;
        type IntoIter = ADIter3<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIter3 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a mut AD3 {
        type Item = f64;
        type IntoIter = ADIterMut3<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIterMut3 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl IntoIterator for AD4 {
        type Item = f64;
        type IntoIter = ADIntoIter4;
        fn into_iter(self) -> Self::IntoIter {
            ADIntoIter4 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a AD4 {
        type Item = f64;
        type IntoIter = ADIter4<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIter4 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a mut AD4 {
        type Item = f64;
        type IntoIter = ADIterMut4<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIterMut4 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl IntoIterator for AD5 {
        type Item = f64;
        type IntoIter = ADIntoIter5;
        fn into_iter(self) -> Self::IntoIter {
            ADIntoIter5 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a AD5 {
        type Item = f64;
        type IntoIter = ADIter5<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIter5 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl<'a> IntoIterator for &'a mut AD5 {
        type Item = f64;
        type IntoIter = ADIterMut5<'a>;
        fn into_iter(self) -> Self::IntoIter {
            ADIterMut5 {
                ad: self,
                index: 0,
                r_index: 0,
            }
        }
    }
    impl Iterator for ADIntoIter1 {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 2 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIter1<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 2 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIterMut1<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 2 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl Iterator for ADIntoIter2 {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 3 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIter2<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 3 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIterMut2<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 3 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl Iterator for ADIntoIter3 {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 4 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIter3<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 4 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIterMut3<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 4 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl Iterator for ADIntoIter4 {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 5 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    4 => self.ad.d4,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIter4<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 5 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    4 => self.ad.d4,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIterMut4<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 5 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    4 => self.ad.d4,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl Iterator for ADIntoIter5 {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 6 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    4 => self.ad.d4,
                    5 => self.ad.d5,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIter5<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 6 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    4 => self.ad.d4,
                    5 => self.ad.d5,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl<'a> Iterator for ADIterMut5<'a> {
        type Item = f64;
        fn next(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index < 6 {
                let result = match self.index {
                    0 => self.ad.d0,
                    1 => self.ad.d1,
                    2 => self.ad.d2,
                    3 => self.ad.d3,
                    4 => self.ad.d4,
                    5 => self.ad.d5,
                    _ => return None,
                };
                self.index += 1;
                Some(result)
            } else {
                None
            }
        }
        fn size_hint(&self) -> (usize, Option<usize>) {
            let lower = self.ad.len() - (self.index + self.r_index);
            let upper = self.ad.len() - (self.index + self.r_index);
            (lower, Some(upper))
        }
    }
    impl FromIterator<f64> for AD1 {
        fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
            let mut z = Self::default();
            for (i, elem) in iter.into_iter().enumerate() {
                z[i] = elem;
            }
            z
        }
    }
    impl FromIterator<f64> for AD2 {
        fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
            let mut z = Self::default();
            for (i, elem) in iter.into_iter().enumerate() {
                z[i] = elem;
            }
            z
        }
    }
    impl FromIterator<f64> for AD3 {
        fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
            let mut z = Self::default();
            for (i, elem) in iter.into_iter().enumerate() {
                z[i] = elem;
            }
            z
        }
    }
    impl FromIterator<f64> for AD4 {
        fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
            let mut z = Self::default();
            for (i, elem) in iter.into_iter().enumerate() {
                z[i] = elem;
            }
            z
        }
    }
    impl FromIterator<f64> for AD5 {
        fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
            let mut z = Self::default();
            for (i, elem) in iter.into_iter().enumerate() {
                z[i] = elem;
            }
            z
        }
    }
    impl DoubleEndedIterator for ADIntoIter1 {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 2 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d1,
                1 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIter1<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 2 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d1,
                1 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIterMut1<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 2 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d1,
                1 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl DoubleEndedIterator for ADIntoIter2 {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 3 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d2,
                1 => self.ad.d1,
                2 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIter2<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 3 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d2,
                1 => self.ad.d1,
                2 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIterMut2<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 3 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d2,
                1 => self.ad.d1,
                2 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl DoubleEndedIterator for ADIntoIter3 {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 4 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d3,
                1 => self.ad.d2,
                2 => self.ad.d1,
                3 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIter3<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 4 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d3,
                1 => self.ad.d2,
                2 => self.ad.d1,
                3 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIterMut3<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 4 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d3,
                1 => self.ad.d2,
                2 => self.ad.d1,
                3 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl DoubleEndedIterator for ADIntoIter4 {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 5 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d4,
                1 => self.ad.d3,
                2 => self.ad.d2,
                3 => self.ad.d1,
                4 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIter4<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 5 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d4,
                1 => self.ad.d3,
                2 => self.ad.d2,
                3 => self.ad.d1,
                4 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIterMut4<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 5 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d4,
                1 => self.ad.d3,
                2 => self.ad.d2,
                3 => self.ad.d1,
                4 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl DoubleEndedIterator for ADIntoIter5 {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 6 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d5,
                1 => self.ad.d4,
                2 => self.ad.d3,
                3 => self.ad.d2,
                4 => self.ad.d1,
                5 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIter5<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 6 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d5,
                1 => self.ad.d4,
                2 => self.ad.d3,
                3 => self.ad.d2,
                4 => self.ad.d1,
                5 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl<'a> DoubleEndedIterator for ADIterMut5<'a> {
        fn next_back(&mut self) -> Option<Self::Item> {
            if self.index + self.r_index == 6 {
                return None;
            }
            let result = match self.r_index {
                0 => self.ad.d5,
                1 => self.ad.d4,
                2 => self.ad.d3,
                3 => self.ad.d2,
                4 => self.ad.d1,
                5 => self.ad.d0,
                _ => return None,
            };
            self.r_index += 1;
            Some(result)
        }
    }
    impl ExactSizeIterator for ADIntoIter1 {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIter1<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIterMut1<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl ExactSizeIterator for ADIntoIter2 {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIter2<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIterMut2<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl ExactSizeIterator for ADIntoIter3 {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIter3<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIterMut3<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl ExactSizeIterator for ADIntoIter4 {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIter4<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIterMut4<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl ExactSizeIterator for ADIntoIter5 {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIter5<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl<'a> ExactSizeIterator for ADIterMut5<'a> {
        fn len(&self) -> usize {
            self.ad.len() - (self.index + self.r_index)
        }
    }
    impl Index<usize> for AD1 {
        type Output = f64;
        fn index(&self, n: usize) -> &Self::Output {
            match n {
                0 => &self.d0,
                1 => &self.d1,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD1"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl IndexMut<usize> for AD1 {
        fn index_mut(&mut self, n: usize) -> &mut Self::Output {
            match n {
                0 => &mut self.d0,
                1 => &mut self.d1,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD1"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl Index<usize> for AD2 {
        type Output = f64;
        fn index(&self, n: usize) -> &Self::Output {
            match n {
                0 => &self.d0,
                1 => &self.d1,
                2 => &self.d2,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD2"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl IndexMut<usize> for AD2 {
        fn index_mut(&mut self, n: usize) -> &mut Self::Output {
            match n {
                0 => &mut self.d0,
                1 => &mut self.d1,
                2 => &mut self.d2,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD2"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl Index<usize> for AD3 {
        type Output = f64;
        fn index(&self, n: usize) -> &Self::Output {
            match n {
                0 => &self.d0,
                1 => &self.d1,
                2 => &self.d2,
                3 => &self.d3,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD3"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl IndexMut<usize> for AD3 {
        fn index_mut(&mut self, n: usize) -> &mut Self::Output {
            match n {
                0 => &mut self.d0,
                1 => &mut self.d1,
                2 => &mut self.d2,
                3 => &mut self.d3,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD3"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl Index<usize> for AD4 {
        type Output = f64;
        fn index(&self, n: usize) -> &Self::Output {
            match n {
                0 => &self.d0,
                1 => &self.d1,
                2 => &self.d2,
                3 => &self.d3,
                4 => &self.d4,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD4"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl IndexMut<usize> for AD4 {
        fn index_mut(&mut self, n: usize) -> &mut Self::Output {
            match n {
                0 => &mut self.d0,
                1 => &mut self.d1,
                2 => &mut self.d2,
                3 => &mut self.d3,
                4 => &mut self.d4,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD4"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl Index<usize> for AD5 {
        type Output = f64;
        fn index(&self, n: usize) -> &Self::Output {
            match n {
                0 => &self.d0,
                1 => &self.d1,
                2 => &self.d2,
                3 => &self.d3,
                4 => &self.d4,
                5 => &self.d5,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD5"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl IndexMut<usize> for AD5 {
        fn index_mut(&mut self, n: usize) -> &mut Self::Output {
            match n {
                0 => &mut self.d0,
                1 => &mut self.d1,
                2 => &mut self.d2,
                3 => &mut self.d3,
                4 => &mut self.d4,
                5 => &mut self.d5,
                _ => ::std::rt::begin_panic_fmt(&::core::fmt::Arguments::new_v1(
                    &["", " exceed order of AD5"],
                    &match (&n,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::Display::fmt,
                        )],
                    },
                )),
            }
        }
    }
    impl Neg for AD1 {
        type Output = Self;
        fn neg(self) -> Self::Output {
            AD1::new(-self.d0, -self.d1)
        }
    }
    impl Neg for AD2 {
        type Output = Self;
        fn neg(self) -> Self::Output {
            AD2::new(-self.d0, -self.d1, -self.d2)
        }
    }
    impl Neg for AD3 {
        type Output = Self;
        fn neg(self) -> Self::Output {
            AD3::new(-self.d0, -self.d1, -self.d2, -self.d3)
        }
    }
    impl Neg for AD4 {
        type Output = Self;
        fn neg(self) -> Self::Output {
            AD4::new(-self.d0, -self.d1, -self.d2, -self.d3, -self.d4)
        }
    }
    impl Neg for AD5 {
        type Output = Self;
        fn neg(self) -> Self::Output {
            AD5::new(-self.d0, -self.d1, -self.d2, -self.d3, -self.d4, -self.d5)
        }
    }
    impl Add<AD1> for AD1 {
        type Output = AD1;
        fn add(self, rhs: AD1) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Add<AD2> for AD1 {
        type Output = AD2;
        fn add(self, rhs: AD2) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Add<AD3> for AD1 {
        type Output = AD3;
        fn add(self, rhs: AD3) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Add<AD4> for AD1 {
        type Output = AD4;
        fn add(self, rhs: AD4) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Add<AD5> for AD1 {
        type Output = AD5;
        fn add(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Add<AD1> for AD2 {
        type Output = AD2;
        fn add(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z
        }
    }
    impl Add<AD2> for AD2 {
        type Output = AD2;
        fn add(self, rhs: AD2) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z
        }
    }
    impl Add<AD3> for AD2 {
        type Output = AD3;
        fn add(self, rhs: AD3) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z
        }
    }
    impl Add<AD4> for AD2 {
        type Output = AD4;
        fn add(self, rhs: AD4) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z
        }
    }
    impl Add<AD5> for AD2 {
        type Output = AD5;
        fn add(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z
        }
    }
    impl Add<AD1> for AD3 {
        type Output = AD3;
        fn add(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z
        }
    }
    impl Add<AD2> for AD3 {
        type Output = AD3;
        fn add(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z.d2 += rhs.d2;
            z
        }
    }
    impl Add<AD3> for AD3 {
        type Output = AD3;
        fn add(self, rhs: AD3) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z
        }
    }
    impl Add<AD4> for AD3 {
        type Output = AD4;
        fn add(self, rhs: AD4) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z
        }
    }
    impl Add<AD5> for AD3 {
        type Output = AD5;
        fn add(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z
        }
    }
    impl Add<AD1> for AD4 {
        type Output = AD4;
        fn add(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z
        }
    }
    impl Add<AD2> for AD4 {
        type Output = AD4;
        fn add(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z.d2 += rhs.d2;
            z
        }
    }
    impl Add<AD3> for AD4 {
        type Output = AD4;
        fn add(self, rhs: AD3) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z.d2 += rhs.d2;
            z.d3 += rhs.d3;
            z
        }
    }
    impl Add<AD4> for AD4 {
        type Output = AD4;
        fn add(self, rhs: AD4) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z.d4 += self.d4;
            z
        }
    }
    impl Add<AD5> for AD4 {
        type Output = AD5;
        fn add(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z.d4 += self.d4;
            z
        }
    }
    impl Add<AD1> for AD5 {
        type Output = AD5;
        fn add(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z
        }
    }
    impl Add<AD2> for AD5 {
        type Output = AD5;
        fn add(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z.d2 += rhs.d2;
            z
        }
    }
    impl Add<AD3> for AD5 {
        type Output = AD5;
        fn add(self, rhs: AD3) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z.d2 += rhs.d2;
            z.d3 += rhs.d3;
            z
        }
    }
    impl Add<AD4> for AD5 {
        type Output = AD5;
        fn add(self, rhs: AD4) -> Self::Output {
            let mut z = self.clone();
            z.d0 += rhs.d0;
            z.d1 += rhs.d1;
            z.d2 += rhs.d2;
            z.d3 += rhs.d3;
            z.d4 += rhs.d4;
            z
        }
    }
    impl Add<AD5> for AD5 {
        type Output = AD5;
        fn add(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z.d4 += self.d4;
            z.d5 += self.d5;
            z
        }
    }
    impl Sub<AD1> for AD1 {
        type Output = AD1;
        fn sub(self, rhs: AD1) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Sub<AD2> for AD1 {
        type Output = AD2;
        fn sub(self, rhs: AD2) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Sub<AD3> for AD1 {
        type Output = AD3;
        fn sub(self, rhs: AD3) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Sub<AD4> for AD1 {
        type Output = AD4;
        fn sub(self, rhs: AD4) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Sub<AD5> for AD1 {
        type Output = AD5;
        fn sub(self, rhs: AD5) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z
        }
    }
    impl Sub<AD1> for AD2 {
        type Output = AD2;
        fn sub(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z
        }
    }
    impl Sub<AD2> for AD2 {
        type Output = AD2;
        fn sub(self, rhs: AD2) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z
        }
    }
    impl Sub<AD3> for AD2 {
        type Output = AD3;
        fn sub(self, rhs: AD3) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z
        }
    }
    impl Sub<AD4> for AD2 {
        type Output = AD4;
        fn sub(self, rhs: AD4) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z
        }
    }
    impl Sub<AD5> for AD2 {
        type Output = AD5;
        fn sub(self, rhs: AD5) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z
        }
    }
    impl Sub<AD1> for AD3 {
        type Output = AD3;
        fn sub(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z
        }
    }
    impl Sub<AD2> for AD3 {
        type Output = AD3;
        fn sub(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z.d2 -= rhs.d2;
            z
        }
    }
    impl Sub<AD3> for AD3 {
        type Output = AD3;
        fn sub(self, rhs: AD3) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z
        }
    }
    impl Sub<AD4> for AD3 {
        type Output = AD4;
        fn sub(self, rhs: AD4) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z
        }
    }
    impl Sub<AD5> for AD3 {
        type Output = AD5;
        fn sub(self, rhs: AD5) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z
        }
    }
    impl Sub<AD1> for AD4 {
        type Output = AD4;
        fn sub(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z
        }
    }
    impl Sub<AD2> for AD4 {
        type Output = AD4;
        fn sub(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z.d2 -= rhs.d2;
            z
        }
    }
    impl Sub<AD3> for AD4 {
        type Output = AD4;
        fn sub(self, rhs: AD3) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z.d2 -= rhs.d2;
            z.d3 -= rhs.d3;
            z
        }
    }
    impl Sub<AD4> for AD4 {
        type Output = AD4;
        fn sub(self, rhs: AD4) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z.d4 += self.d4;
            z
        }
    }
    impl Sub<AD5> for AD4 {
        type Output = AD5;
        fn sub(self, rhs: AD5) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z.d4 += self.d4;
            z
        }
    }
    impl Sub<AD1> for AD5 {
        type Output = AD5;
        fn sub(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z
        }
    }
    impl Sub<AD2> for AD5 {
        type Output = AD5;
        fn sub(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z.d2 -= rhs.d2;
            z
        }
    }
    impl Sub<AD3> for AD5 {
        type Output = AD5;
        fn sub(self, rhs: AD3) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z.d2 -= rhs.d2;
            z.d3 -= rhs.d3;
            z
        }
    }
    impl Sub<AD4> for AD5 {
        type Output = AD5;
        fn sub(self, rhs: AD4) -> Self::Output {
            let mut z = self.clone();
            z.d0 -= rhs.d0;
            z.d1 -= rhs.d1;
            z.d2 -= rhs.d2;
            z.d3 -= rhs.d3;
            z.d4 -= rhs.d4;
            z
        }
    }
    impl Sub<AD5> for AD5 {
        type Output = AD5;
        fn sub(self, rhs: AD5) -> Self::Output {
            let mut z = -rhs.clone();
            z.d0 += self.d0;
            z.d1 += self.d1;
            z.d2 += self.d2;
            z.d3 += self.d3;
            z.d4 += self.d4;
            z.d5 += self.d5;
            z
        }
    }
    impl Mul<AD1> for AD1 {
        type Output = AD1;
        fn mul(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD1> for AD1 {
        type Output = AD1;
        fn mul(self, rhs: &'a AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD2> for AD1 {
        type Output = AD2;
        fn mul(self, rhs: AD2) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD2> for AD1 {
        type Output = AD2;
        fn mul(self, rhs: &'a AD2) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD3> for AD1 {
        type Output = AD3;
        fn mul(self, rhs: AD3) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD3> for AD1 {
        type Output = AD3;
        fn mul(self, rhs: &'a AD3) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD4> for AD1 {
        type Output = AD4;
        fn mul(self, rhs: AD4) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD4> for AD1 {
        type Output = AD4;
        fn mul(self, rhs: &'a AD4) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD5> for AD1 {
        type Output = AD5;
        fn mul(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD5> for AD1 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD5) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD1> for AD2 {
        type Output = AD2;
        fn mul(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD1> for AD2 {
        type Output = AD2;
        fn mul(self, rhs: &'a AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD2> for AD2 {
        type Output = AD2;
        fn mul(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD2> for AD2 {
        type Output = AD2;
        fn mul(self, rhs: &'a AD2) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD3> for AD2 {
        type Output = AD3;
        fn mul(self, rhs: AD3) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD3> for AD2 {
        type Output = AD3;
        fn mul(self, rhs: &'a AD3) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD4> for AD2 {
        type Output = AD4;
        fn mul(self, rhs: AD4) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD4> for AD2 {
        type Output = AD4;
        fn mul(self, rhs: &'a AD4) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD5> for AD2 {
        type Output = AD5;
        fn mul(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD5> for AD2 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD5) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD1> for AD3 {
        type Output = AD3;
        fn mul(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD1> for AD3 {
        type Output = AD3;
        fn mul(self, rhs: &'a AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD2> for AD3 {
        type Output = AD3;
        fn mul(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD2> for AD3 {
        type Output = AD3;
        fn mul(self, rhs: &'a AD2) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD3> for AD3 {
        type Output = AD3;
        fn mul(self, rhs: AD3) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD3> for AD3 {
        type Output = AD3;
        fn mul(self, rhs: &'a AD3) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD4> for AD3 {
        type Output = AD4;
        fn mul(self, rhs: AD4) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD4> for AD3 {
        type Output = AD4;
        fn mul(self, rhs: &'a AD4) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD5> for AD3 {
        type Output = AD5;
        fn mul(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD5> for AD3 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD5) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD1> for AD4 {
        type Output = AD4;
        fn mul(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD1> for AD4 {
        type Output = AD4;
        fn mul(self, rhs: &'a AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD2> for AD4 {
        type Output = AD4;
        fn mul(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD2> for AD4 {
        type Output = AD4;
        fn mul(self, rhs: &'a AD2) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD3> for AD4 {
        type Output = AD4;
        fn mul(self, rhs: AD3) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD3> for AD4 {
        type Output = AD4;
        fn mul(self, rhs: &'a AD3) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD4> for AD4 {
        type Output = AD4;
        fn mul(self, rhs: AD4) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD4> for AD4 {
        type Output = AD4;
        fn mul(self, rhs: &'a AD4) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD5> for AD4 {
        type Output = AD5;
        fn mul(self, rhs: AD5) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(rhs.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD5> for AD4 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD5) -> Self::Output {
            let mut z = rhs.clone();
            let x = Self::Output::from(self);
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = x
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD1> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD1> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD1) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD2> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: AD2) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD2> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD2) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD3> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: AD3) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD3> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD3) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD4> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: AD4) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD4> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD4) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Mul<AD5> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: AD5) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl<'a> Mul<&'a AD5> for AD5 {
        type Output = AD5;
        fn mul(self, rhs: &'a AD5) -> Self::Output {
            let mut z = self.clone();
            let y = Self::Output::from(rhs);
            for t in 0..z.len() {
                z[t] = self
                    .iter()
                    .take(t + 1)
                    .zip(y.iter().take(t + 1).rev())
                    .enumerate()
                    .fold(0f64, |s, (k, (x1, y1))| s + (C(t, k) as f64) * x1 * y1)
            }
            z
        }
    }
    impl Div<AD1> for AD1 {
        type Output = AD1;
        fn div(self, rhs: AD1) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD2> for AD1 {
        type Output = AD2;
        fn div(self, rhs: AD2) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD3> for AD1 {
        type Output = AD3;
        fn div(self, rhs: AD3) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD4> for AD1 {
        type Output = AD4;
        fn div(self, rhs: AD4) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD5> for AD1 {
        type Output = AD5;
        fn div(self, rhs: AD5) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD1> for AD2 {
        type Output = AD2;
        fn div(self, rhs: AD1) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD2> for AD2 {
        type Output = AD2;
        fn div(self, rhs: AD2) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD3> for AD2 {
        type Output = AD3;
        fn div(self, rhs: AD3) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD4> for AD2 {
        type Output = AD4;
        fn div(self, rhs: AD4) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD5> for AD2 {
        type Output = AD5;
        fn div(self, rhs: AD5) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD1> for AD3 {
        type Output = AD3;
        fn div(self, rhs: AD1) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD2> for AD3 {
        type Output = AD3;
        fn div(self, rhs: AD2) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD3> for AD3 {
        type Output = AD3;
        fn div(self, rhs: AD3) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD4> for AD3 {
        type Output = AD4;
        fn div(self, rhs: AD4) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD5> for AD3 {
        type Output = AD5;
        fn div(self, rhs: AD5) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD1> for AD4 {
        type Output = AD4;
        fn div(self, rhs: AD1) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD2> for AD4 {
        type Output = AD4;
        fn div(self, rhs: AD2) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD3> for AD4 {
        type Output = AD4;
        fn div(self, rhs: AD3) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD4> for AD4 {
        type Output = AD4;
        fn div(self, rhs: AD4) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD5> for AD4 {
        type Output = AD5;
        fn div(self, rhs: AD5) -> Self::Output {
            let mut z = Self::Output::default();
            let x = Self::Output::from(self);
            z[0] = x[0] / rhs[0];
            let y0 = 1f64 / rhs[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in rhs
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (x[i] - s);
            }
            z
        }
    }
    impl Div<AD1> for AD5 {
        type Output = AD5;
        fn div(self, rhs: AD1) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD2> for AD5 {
        type Output = AD5;
        fn div(self, rhs: AD2) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD3> for AD5 {
        type Output = AD5;
        fn div(self, rhs: AD3) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD4> for AD5 {
        type Output = AD5;
        fn div(self, rhs: AD4) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl Div<AD5> for AD5 {
        type Output = AD5;
        fn div(self, rhs: AD5) -> Self::Output {
            let mut z = Self::Output::default();
            let y = Self::Output::from(rhs);
            z[0] = self[0] / y[0];
            let y0 = 1f64 / y[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (y1, z1)) in y
                    .iter()
                    .skip(1)
                    .take(i)
                    .zip(z.iter().take(i).rev())
                    .enumerate()
                {
                    s += (C(i, j + 1) as f64) * y1 * z1;
                }
                z[i] = y0 * (self[i] - s);
            }
            z
        }
    }
    impl ExpLogOps for AD1 {
        fn exp(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].exp();
            for i in 1..z.len() {
                z[i] = z
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (z1, x1))| x + (C(i - 1, k) as f64) * x1 * z1);
            }
            z
        }
        fn ln(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].ln();
            let x0 = 1f64 / self[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (k, (z1, x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(self.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, k + 1) as f64) * z1 * x1;
                }
                z[i] = x0 * (self[i] - s);
            }
            z
        }
        fn log(&self, base: f64) -> Self {
            self.ln().iter().map(|x| x / base.ln()).collect()
        }
    }
    impl ExpLogOps for AD2 {
        fn exp(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].exp();
            for i in 1..z.len() {
                z[i] = z
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (z1, x1))| x + (C(i - 1, k) as f64) * x1 * z1);
            }
            z
        }
        fn ln(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].ln();
            let x0 = 1f64 / self[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (k, (z1, x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(self.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, k + 1) as f64) * z1 * x1;
                }
                z[i] = x0 * (self[i] - s);
            }
            z
        }
        fn log(&self, base: f64) -> Self {
            self.ln().iter().map(|x| x / base.ln()).collect()
        }
    }
    impl ExpLogOps for AD3 {
        fn exp(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].exp();
            for i in 1..z.len() {
                z[i] = z
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (z1, x1))| x + (C(i - 1, k) as f64) * x1 * z1);
            }
            z
        }
        fn ln(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].ln();
            let x0 = 1f64 / self[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (k, (z1, x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(self.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, k + 1) as f64) * z1 * x1;
                }
                z[i] = x0 * (self[i] - s);
            }
            z
        }
        fn log(&self, base: f64) -> Self {
            self.ln().iter().map(|x| x / base.ln()).collect()
        }
    }
    impl ExpLogOps for AD4 {
        fn exp(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].exp();
            for i in 1..z.len() {
                z[i] = z
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (z1, x1))| x + (C(i - 1, k) as f64) * x1 * z1);
            }
            z
        }
        fn ln(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].ln();
            let x0 = 1f64 / self[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (k, (z1, x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(self.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, k + 1) as f64) * z1 * x1;
                }
                z[i] = x0 * (self[i] - s);
            }
            z
        }
        fn log(&self, base: f64) -> Self {
            self.ln().iter().map(|x| x / base.ln()).collect()
        }
    }
    impl ExpLogOps for AD5 {
        fn exp(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].exp();
            for i in 1..z.len() {
                z[i] = z
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (z1, x1))| x + (C(i - 1, k) as f64) * x1 * z1);
            }
            z
        }
        fn ln(&self) -> Self {
            let mut z = Self::default();
            z[0] = self[0].ln();
            let x0 = 1f64 / self[0];
            for i in 1..z.len() {
                let mut s = 0f64;
                for (k, (z1, x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(self.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, k + 1) as f64) * z1 * x1;
                }
                z[i] = x0 * (self[i] - s);
            }
            z
        }
        fn log(&self, base: f64) -> Self {
            self.ln().iter().map(|x| x / base.ln()).collect()
        }
    }
    impl PowOps for AD1 {
        fn powi(&self, n: i32) -> Self {
            let mut z = self.clone();
            for _i in 1..n {
                z = z * self;
            }
            z
        }
        fn powf(&self, f: f64) -> Self {
            let ln_x = self.ln();
            let mut z = Self::default();
            z[0] = self.d0.powf(f);
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (z1, ln_x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(ln_x.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, j + 1) as f64) * z1 * ln_x1;
                }
                z[i] = f * (z[0] * ln_x[i] + s);
            }
            z
        }
        fn pow(&self, _f: Self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl PowOps for AD2 {
        fn powi(&self, n: i32) -> Self {
            let mut z = self.clone();
            for _i in 1..n {
                z = z * self;
            }
            z
        }
        fn powf(&self, f: f64) -> Self {
            let ln_x = self.ln();
            let mut z = Self::default();
            z[0] = self.d0.powf(f);
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (z1, ln_x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(ln_x.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, j + 1) as f64) * z1 * ln_x1;
                }
                z[i] = f * (z[0] * ln_x[i] + s);
            }
            z
        }
        fn pow(&self, _f: Self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl PowOps for AD3 {
        fn powi(&self, n: i32) -> Self {
            let mut z = self.clone();
            for _i in 1..n {
                z = z * self;
            }
            z
        }
        fn powf(&self, f: f64) -> Self {
            let ln_x = self.ln();
            let mut z = Self::default();
            z[0] = self.d0.powf(f);
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (z1, ln_x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(ln_x.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, j + 1) as f64) * z1 * ln_x1;
                }
                z[i] = f * (z[0] * ln_x[i] + s);
            }
            z
        }
        fn pow(&self, _f: Self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl PowOps for AD4 {
        fn powi(&self, n: i32) -> Self {
            let mut z = self.clone();
            for _i in 1..n {
                z = z * self;
            }
            z
        }
        fn powf(&self, f: f64) -> Self {
            let ln_x = self.ln();
            let mut z = Self::default();
            z[0] = self.d0.powf(f);
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (z1, ln_x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(ln_x.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, j + 1) as f64) * z1 * ln_x1;
                }
                z[i] = f * (z[0] * ln_x[i] + s);
            }
            z
        }
        fn pow(&self, _f: Self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl PowOps for AD5 {
        fn powi(&self, n: i32) -> Self {
            let mut z = self.clone();
            for _i in 1..n {
                z = z * self;
            }
            z
        }
        fn powf(&self, f: f64) -> Self {
            let ln_x = self.ln();
            let mut z = Self::default();
            z[0] = self.d0.powf(f);
            for i in 1..z.len() {
                let mut s = 0f64;
                for (j, (z1, ln_x1)) in z
                    .iter()
                    .skip(1)
                    .take(i - 1)
                    .zip(ln_x.iter().skip(1).take(i - 1).rev())
                    .enumerate()
                {
                    s += (C(i - 1, j + 1) as f64) * z1 * ln_x1;
                }
                z[i] = f * (z[0] * ln_x[i] + s);
            }
            z
        }
        fn pow(&self, _f: Self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl TrigOps for AD1 {
        fn sin_cos(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sin();
            v[0] = self[0].cos();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = -u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn sinh_cosh(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sinh();
            v[0] = self[0].cosh();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn asin(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acos(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atan(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn asinh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acosh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atanh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl TrigOps for AD2 {
        fn sin_cos(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sin();
            v[0] = self[0].cos();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = -u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn sinh_cosh(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sinh();
            v[0] = self[0].cosh();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn asin(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acos(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atan(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn asinh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acosh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atanh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl TrigOps for AD3 {
        fn sin_cos(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sin();
            v[0] = self[0].cos();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = -u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn sinh_cosh(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sinh();
            v[0] = self[0].cosh();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn asin(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acos(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atan(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn asinh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acosh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atanh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl TrigOps for AD4 {
        fn sin_cos(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sin();
            v[0] = self[0].cos();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = -u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn sinh_cosh(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sinh();
            v[0] = self[0].cosh();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn asin(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acos(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atan(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn asinh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acosh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atanh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl TrigOps for AD5 {
        fn sin_cos(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sin();
            v[0] = self[0].cos();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = -u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn sinh_cosh(&self) -> (Self, Self) {
            let mut u = Self::default();
            let mut v = Self::default();
            u[0] = self[0].sinh();
            v[0] = self[0].cosh();
            for i in 1..u.len() {
                u[i] = v
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (v1, x1))| x + (C(i - 1, k) as f64) * x1 * v1);
                v[i] = u
                    .iter()
                    .take(i)
                    .zip(self.iter().skip(1).take(i).rev())
                    .enumerate()
                    .fold(0f64, |x, (k, (u1, x1))| x + (C(i - 1, k) as f64) * x1 * u1);
            }
            (u, v)
        }
        fn asin(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acos(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atan(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn asinh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn acosh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
        fn atanh(&self) -> Self {
            {
                ::std::rt::begin_panic("not implemented")
            }
        }
    }
    impl AD for AD1 {}
    impl AD for AD2 {}
    impl AD for AD3 {}
    impl AD for AD4 {}
    impl AD for AD5 {}
    pub trait AD:
        std::fmt::Display
        + Clone
        + Copy
        + PartialEq
        + From<f64>
        + Into<f64>
        + Add<Output = Self>
        + Sub<Output = Self>
        + Mul<Output = Self>
        + Div<Output = Self>
        + Add<f64, Output = Self>
        + Sub<f64, Output = Self>
        + Mul<f64, Output = Self>
        + Div<f64, Output = Self>
        + PowOps
        + ExpLogOps
        + TrigOps
        + From<AD1>
        + Into<AD1>
        + From<AD2>
        + Into<AD2>
        + From<AD3>
        + Into<AD3>
        + From<AD4>
        + Into<AD4>
        + From<AD5>
        + Into<AD5>
    {
        fn to_ad1(self) -> AD1 {
            self.into()
        }
        fn to_ad2(self) -> AD2 {
            self.into()
        }
        fn to_ad3(self) -> AD3 {
            self.into()
        }
        fn to_ad4(self) -> AD4 {
            self.into()
        }
        fn to_ad5(self) -> AD5 {
            self.into()
        }
    }
    impl From<f64> for AD1 {
        fn from(other: f64) -> Self {
            let f = other as f64;
            let mut z = Self::default();
            z.d0 = f;
            z
        }
    }
    impl From<f64> for AD2 {
        fn from(other: f64) -> Self {
            let f = other as f64;
            let mut z = Self::default();
            z.d0 = f;
            z
        }
    }
    impl From<f64> for AD3 {
        fn from(other: f64) -> Self {
            let f = other as f64;
            let mut z = Self::default();
            z.d0 = f;
            z
        }
    }
    impl From<f64> for AD4 {
        fn from(other: f64) -> Self {
            let f = other as f64;
            let mut z = Self::default();
            z.d0 = f;
            z
        }
    }
    impl From<f64> for AD5 {
        fn from(other: f64) -> Self {
            let f = other as f64;
            let mut z = Self::default();
            z.d0 = f;
            z
        }
    }
    impl Add<f64> for AD1 {
        type Output = Self;
        fn add(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 += rhs;
            z
        }
    }
    impl Add<f64> for AD2 {
        type Output = Self;
        fn add(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 += rhs;
            z
        }
    }
    impl Add<f64> for AD3 {
        type Output = Self;
        fn add(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 += rhs;
            z
        }
    }
    impl Add<f64> for AD4 {
        type Output = Self;
        fn add(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 += rhs;
            z
        }
    }
    impl Add<f64> for AD5 {
        type Output = Self;
        fn add(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 += rhs;
            z
        }
    }
    impl Sub<f64> for AD1 {
        type Output = Self;
        fn sub(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 -= rhs;
            z
        }
    }
    impl Sub<f64> for AD2 {
        type Output = Self;
        fn sub(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 -= rhs;
            z
        }
    }
    impl Sub<f64> for AD3 {
        type Output = Self;
        fn sub(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 -= rhs;
            z
        }
    }
    impl Sub<f64> for AD4 {
        type Output = Self;
        fn sub(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 -= rhs;
            z
        }
    }
    impl Sub<f64> for AD5 {
        type Output = Self;
        fn sub(self, rhs: f64) -> Self::Output {
            let mut z = self;
            z.d0 -= rhs;
            z
        }
    }
    impl Mul<f64> for AD1 {
        type Output = Self;
        fn mul(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x * rhs).collect()
        }
    }
    impl Mul<f64> for AD2 {
        type Output = Self;
        fn mul(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x * rhs).collect()
        }
    }
    impl Mul<f64> for AD3 {
        type Output = Self;
        fn mul(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x * rhs).collect()
        }
    }
    impl Mul<f64> for AD4 {
        type Output = Self;
        fn mul(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x * rhs).collect()
        }
    }
    impl Mul<f64> for AD5 {
        type Output = Self;
        fn mul(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x * rhs).collect()
        }
    }
    impl Div<f64> for AD1 {
        type Output = Self;
        fn div(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x / rhs).collect()
        }
    }
    impl Div<f64> for AD2 {
        type Output = Self;
        fn div(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x / rhs).collect()
        }
    }
    impl Div<f64> for AD3 {
        type Output = Self;
        fn div(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x / rhs).collect()
        }
    }
    impl Div<f64> for AD4 {
        type Output = Self;
        fn div(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x / rhs).collect()
        }
    }
    impl Div<f64> for AD5 {
        type Output = Self;
        fn div(self, rhs: f64) -> Self::Output {
            self.iter().map(|x| x / rhs).collect()
        }
    }
    impl Add<AD1> for f64 {
        type Output = AD1;
        fn add(self, rhs: AD1) -> Self::Output {
            let mut z = rhs;
            z.d0 += self;
            z
        }
    }
    impl Add<AD2> for f64 {
        type Output = AD2;
        fn add(self, rhs: AD2) -> Self::Output {
            let mut z = rhs;
            z.d0 += self;
            z
        }
    }
    impl Add<AD3> for f64 {
        type Output = AD3;
        fn add(self, rhs: AD3) -> Self::Output {
            let mut z = rhs;
            z.d0 += self;
            z
        }
    }
    impl Add<AD4> for f64 {
        type Output = AD4;
        fn add(self, rhs: AD4) -> Self::Output {
            let mut z = rhs;
            z.d0 += self;
            z
        }
    }
    impl Add<AD5> for f64 {
        type Output = AD5;
        fn add(self, rhs: AD5) -> Self::Output {
            let mut z = rhs;
            z.d0 += self;
            z
        }
    }
    impl Sub<AD1> for f64 {
        type Output = AD1;
        fn sub(self, rhs: AD1) -> Self::Output {
            let mut z = -rhs;
            z.d0 += self;
            z
        }
    }
    impl Sub<AD2> for f64 {
        type Output = AD2;
        fn sub(self, rhs: AD2) -> Self::Output {
            let mut z = -rhs;
            z.d0 += self;
            z
        }
    }
    impl Sub<AD3> for f64 {
        type Output = AD3;
        fn sub(self, rhs: AD3) -> Self::Output {
            let mut z = -rhs;
            z.d0 += self;
            z
        }
    }
    impl Sub<AD4> for f64 {
        type Output = AD4;
        fn sub(self, rhs: AD4) -> Self::Output {
            let mut z = -rhs;
            z.d0 += self;
            z
        }
    }
    impl Sub<AD5> for f64 {
        type Output = AD5;
        fn sub(self, rhs: AD5) -> Self::Output {
            let mut z = -rhs;
            z.d0 += self;
            z
        }
    }
    impl Mul<AD1> for f64 {
        type Output = AD1;
        fn mul(self, rhs: AD1) -> Self::Output {
            rhs.iter().map(|x| x * self).collect()
        }
    }
    impl Mul<AD2> for f64 {
        type Output = AD2;
        fn mul(self, rhs: AD2) -> Self::Output {
            rhs.iter().map(|x| x * self).collect()
        }
    }
    impl Mul<AD3> for f64 {
        type Output = AD3;
        fn mul(self, rhs: AD3) -> Self::Output {
            rhs.iter().map(|x| x * self).collect()
        }
    }
    impl Mul<AD4> for f64 {
        type Output = AD4;
        fn mul(self, rhs: AD4) -> Self::Output {
            rhs.iter().map(|x| x * self).collect()
        }
    }
    impl Mul<AD5> for f64 {
        type Output = AD5;
        fn mul(self, rhs: AD5) -> Self::Output {
            rhs.iter().map(|x| x * self).collect()
        }
    }
    impl Div<AD1> for f64 {
        type Output = AD1;
        fn div(self, rhs: AD1) -> Self::Output {
            let ad1 = AD1::from(self);
            ad1 / rhs
        }
    }
    impl Div<AD2> for f64 {
        type Output = AD2;
        fn div(self, rhs: AD2) -> Self::Output {
            let ad1 = AD1::from(self);
            ad1 / rhs
        }
    }
    impl Div<AD3> for f64 {
        type Output = AD3;
        fn div(self, rhs: AD3) -> Self::Output {
            let ad1 = AD1::from(self);
            ad1 / rhs
        }
    }
    impl Div<AD4> for f64 {
        type Output = AD4;
        fn div(self, rhs: AD4) -> Self::Output {
            let ad1 = AD1::from(self);
            ad1 / rhs
        }
    }
    impl Div<AD5> for f64 {
        type Output = AD5;
        fn div(self, rhs: AD5) -> Self::Output {
            let ad1 = AD1::from(self);
            ad1 / rhs
        }
    }
    impl From<AD1> for f64 {
        fn from(ad: AD1) -> Self {
            ad.d0
        }
    }
    impl From<AD2> for f64 {
        fn from(ad: AD2) -> Self {
            ad.d0
        }
    }
    impl From<AD3> for f64 {
        fn from(ad: AD3) -> Self {
            ad.d0
        }
    }
    impl From<AD4> for f64 {
        fn from(ad: AD4) -> Self {
            ad.d0
        }
    }
    impl From<AD5> for f64 {
        fn from(ad: AD5) -> Self {
            ad.d0
        }
    }
    impl<F: Fn(AD2) -> AD2> StableFn<AD2> for ADLift<F, AD2> {
        type Output = AD2;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD2::from(target)).into()
        }
    }
    impl<F: Fn(AD3) -> AD3> StableFn<AD2> for ADLift<F, AD3> {
        type Output = AD2;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD3::from(target)).into()
        }
    }
    impl<F: Fn(AD3) -> AD3> StableFn<AD3> for ADLift<F, AD3> {
        type Output = AD3;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD3::from(target)).into()
        }
    }
    impl<F: Fn(AD4) -> AD4> StableFn<AD2> for ADLift<F, AD4> {
        type Output = AD2;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD4::from(target)).into()
        }
    }
    impl<F: Fn(AD4) -> AD4> StableFn<AD3> for ADLift<F, AD4> {
        type Output = AD3;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD4::from(target)).into()
        }
    }
    impl<F: Fn(AD4) -> AD4> StableFn<AD4> for ADLift<F, AD4> {
        type Output = AD4;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD4::from(target)).into()
        }
    }
    impl<F: Fn(AD5) -> AD5> StableFn<AD2> for ADLift<F, AD5> {
        type Output = AD2;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD5::from(target)).into()
        }
    }
    impl<F: Fn(AD5) -> AD5> StableFn<AD3> for ADLift<F, AD5> {
        type Output = AD3;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD5::from(target)).into()
        }
    }
    impl<F: Fn(AD5) -> AD5> StableFn<AD4> for ADLift<F, AD5> {
        type Output = AD4;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD5::from(target)).into()
        }
    }
    impl<F: Fn(AD5) -> AD5> StableFn<AD5> for ADLift<F, AD5> {
        type Output = AD5;
        fn call_stable(&self, target: Self::Output) -> Self::Output {
            self.f(AD5::from(target)).into()
        }
    }
    /// Lift AD functions
    ///
    /// # Description
    /// To lift `AD` functions
    ///
    /// # Implementation
    ///
    /// * All `Fn(T) -> T where T:AD` functions can be lift to `Fn(f64) -> f64`
    /// * If `j > i`, then `Fn(AD{j}) -> AD{j}` can be lift to `Fn(AD{i}) -> AD{i}`
    ///
    /// # Usage
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let ad0 = 2f64;
    ///     let ad1 = AD1::new(2f64, 1f64);
    ///     
    ///     let lift1 = ADLift::<_, AD1>::new(f_ad);
    ///     let lift2 = ADLift::<_, AD2>::new(f_ad2);
    ///
    ///     let ans_ad0 = ad0.powi(2);
    ///     let ans_ad1 = ad1.powi(2);
    ///
    ///     // All AD function can be lift to f64
    ///     assert_eq!(ans_ad0, lift1.call_stable(ad0));
    ///     assert_eq!(ans_ad0, lift2.call_stable(ad0));
    ///
    ///     // AD2 is higher than AD1 (AD2 -> AD1 lifting is allowed)
    ///     assert_eq!(ans_ad1, lift2.call_stable(ad1));
    /// }
    ///
    /// fn f_ad<T: AD>(x: T) -> T {
    ///     x.powi(2)
    /// }
    ///
    /// fn f_ad2(x: AD2) -> AD2 {
    ///     x.powi(2)
    /// }
    /// ```
    pub struct ADLift<F, T> {
        f: Box<F>,
        _marker: PhantomData<T>,
    }
    impl<F: Fn(T) -> T, T> ADLift<F, T> {
        pub fn new(f: F) -> Self {
            Self {
                f: Box::new(f),
                _marker: PhantomData,
            }
        }
        pub fn f(&self, t: T) -> T {
            (self.f)(t)
        }
    }
    impl<F: Fn(T) -> T, T: AD> StableFn<f64> for ADLift<F, T> {
        type Output = f64;
        fn call_stable(&self, target: f64) -> Self::Output {
            self.f(T::from(target)).into()
        }
    }
    impl<F: Fn(T) -> T, T: AD> StableFn<AD1> for ADLift<F, T> {
        type Output = AD1;
        fn call_stable(&self, target: AD1) -> Self::Output {
            self.f(T::from(target)).into()
        }
    }
}
