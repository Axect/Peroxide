use crate::structure::ad::{AD, AD1, AD2};
use RootState::{P, I};
use std::fmt::Display;
//use std::collections::HashMap;
//use std::marker::PhantomData;

#[derive(Debug, Copy, Clone)]
pub enum RootState {
    P(f64),
    I(f64, f64),
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum RootFind {
    Bisection,
    Newton,
    Secant,
    FalsePosition,
}

/// Structure for Root finding
///
/// * **Caution**: `T` $\geq$ `AD1`
#[derive(Debug, Clone)]
pub struct RootFinder<T: AD> {
    init: RootState,
    pub curr: RootState, 
    method: RootFind,
    f: fn(T) -> T,
    find: bool,
}

#[derive(Debug, Copy, Clone)]
pub enum RootError {
    MismatchedState,
}

impl Display for RootError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            RootError::MismatchedState => {
                writeln!(f, "Mismatched RootState with RootFind method")
            }
        }
    }
}

impl std::error::Error for RootError {}

impl<T: AD> RootFinder<T> {
    pub fn new(init: RootState, method: RootFind, f: fn(T) -> T) -> Result<Self, RootError> {
        match method {
            RootFind::Bisection => {
                match init {
                    P(_) => Err(RootError::MismatchedState),
                    _ => Ok(
                        RootFinder {
                            init,
                            curr: init,
                            method,
                            f,
                            find: false
                        }
                    )
                }
            },
            RootFind::Newton => {
                match init {
                    P(_) => Ok(
                        RootFinder {
                            init,
                            curr: init,
                            method,
                            f,
                            find: false,
                        }
                    ),
                    _ => Err(RootError::MismatchedState),
                }
            },
            RootFind::Secant => {
                match init {
                    P(_) => Err(RootError::MismatchedState),
                    _ => Ok(
                        RootFinder {
                            init,
                            curr: init,
                            method,
                            f,
                            find: false,
                        }
                    )
                }
            },
            RootFind::FalsePosition => {
                match init {
                    P(_) => Ok(
                        RootFinder {
                            init,
                            curr: init,
                            method,
                            f,
                            find: false,
                        }
                    ),
                    _ => Err(RootError::MismatchedState),
                }
            }
        }

    }

    pub fn condition_number(&self) -> f64 {
        match self.curr {
            P(p) => {
                let mut z = AD1::from(p);
                z.d1 = 1f64;
                let fz = (self.f)(z.into()).to_ad1();
                p * fz.d1 / fz.d0
            }
            I(a, b) => {
                let p = (a + b) / 2f64;
                let mut z = AD1::from(p);
                z.d1 = 1f64;
                let fz = (self.f)(z.into()).to_ad1();
                p * fz.d1 / fz.d0
            }
        }
    }

    //#[inline]
    //pub fn update(&mut self) {
    //    match self.method {
    //        RootFind::Bisection => {
    //            match self.curr {
    //                I(a, b) => {
    //                    let x = 0.5 * (a + b);
    //                    if (self.f)(a) * (self.f)(x) < 0f64 {
    //                        self.curr = I(a, x);
    //                    } else if (self.f)(x) * (self.f)(b) < 0f64 {
    //                        self.curr = I(x, b);
    //                    } else {
    //                        self.find = true;
    //                    }
    //                }
    //                _ => unreachable!()
    //            }    
    //        }
    //    }
    //}
}

//pub trait RootFinder {
//    type InitialState;
//    type RootType;
//
//    fn mut_update(&mut self);
//    fn find_root(&mut self) -> Option<Self::RootType>;
//    fn set_initial_condition(&mut self, init: Self::InitialState) -> &mut Self;
//    fn set_times(&mut self) -> &mut Self;
//    fn check_enough(&self) -> bool;
//}
//
//#[derive(Debug, Clone, Copy, Hash, PartialOrd, PartialEq, Eq)]
//enum RootOptions {
//    InitCondition,
//    Times,
//}
//
//pub struct Bisection<F: Fn(f64) -> f64> {
//    init: (f64, f64),
//    curr: (f64, f64),
//    func: Box<F>,
//    times: usize,
//    options: HashMap<RootOptions, bool>,
//}
