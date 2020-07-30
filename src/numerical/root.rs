use crate::traits::stable::StableFn;
use crate::structure::ad::{AD, AD1, AD2, ADLift};
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
    find: RootBool,
    times: usize,
    tol: f64,
    root: f64,
}

#[derive(Debug, Copy, Clone)]
pub enum RootError {
    MismatchedState,
    TimesUp,
    NaNRoot,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum RootBool {
    Find,
    NotYet,
    Error,
}

impl Display for RootError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            RootError::MismatchedState => {
                writeln!(f, "Mismatched RootState with RootFind method")
            },
            RootError::TimesUp => {
                writeln!(f, "No root until set times - More times required")
            },
            RootError::NaNRoot => {
                writeln!(f, "Root is NaN - Should modify initial states")
            }
        }
    }
}

impl std::error::Error for RootError {}

impl<T: AD> RootFinder<T> {
    /// Create RootFinder
    ///
    /// # Usage
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let init = RootState::I(1f64, 3f64);
    ///     let a = RootFinder::new(init, Bisection, f)?;
    ///     println!("{:?}", a);
    /// }
    ///
    /// fn f<T: AD>(x: T) -> T {
    ///     x.powi(2) - 4f64
    /// }
    /// ```
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
                            find: RootBool::NotYet,
                            times: 100,
                            tol: 1e-15,
                            root: 0f64,
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
                            find: RootBool::NotYet,
                            times: 100,
                            tol: 1e-15,
                            root: 0f64,
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
                            find: RootBool::NotYet,
                            times: 100,
                            tol: 1e-15,
                            root: 0f64,
                        }
                    )
                }
            },
            RootFind::FalsePosition => {
                match init {
                    I(_, _) => Ok(
                        RootFinder {
                            init,
                            curr: init,
                            method,
                            f,
                            find: RootBool::NotYet,
                            times: 100,
                            tol: 1e-15,
                            root: 0f64,
                        }
                    ),
                    _ => Err(RootError::MismatchedState),
                }
            }
        }
    }

    /// Condition number
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

    pub fn set_times(&mut self, times: usize) -> &mut Self {
        self.times = times;
        self
    }

    pub fn set_tol(&mut self, tol: f64) -> &mut Self {
        self.tol = tol;
        self
    }

    #[inline]
    pub fn update(&mut self) {
        let lift = ADLift::new(self.f);
        match self.method {
            RootFind::Bisection => {
                match self.curr {
                    I(a, b) => {
                        let x = 0.5 * (a + b);
                        let fa = lift.call_stable(a);
                        let fx = lift.call_stable(x);
                        let fb = lift.call_stable(b);
                        if (a - b).abs() <= self.tol {
                            self.find = RootBool::Find;
                            self.root = x;
                        } else {
                            if fa * fx < 0f64 {
                                self.curr = I(a, x);
                            } else if fx * fb < 0f64 {
                                self.curr = I(x, b);
                            } else if fx == 0f64 {
                                self.find = RootBool::Find;
                                self.root = x;
                            } else {
                                
                            }
                        }
                    }
                    _ => unreachable!()
                }    
            },
            RootFind::FalsePosition => {
                match self.curr {
                    I(a, b) => {
                        let fa = lift.call_stable(a);
                        let fb = lift.call_stable(b);
                        let x = (a * fb - b * fa) / (fb - fa);
                        let fx = lift.call_stable(x);
                        if (a - b).abs() <= self.tol || fx.abs() <= self.tol {
                            self.find = RootBool::Find;
                            self.root = x;
                        } else {
                            if fx * fa > 0f64 {
                                self.curr = I(x, b);
                            } else if fx * fb > 0f64  {
                                self.curr = I(a, x);
                            } else {
                                self.find = RootBool::Find;
                                self.root = x;
                            }
                        }
                    }
                    _ => unreachable!()
                }
            },
            RootFind::Newton => {
                unimplemented!()
            },
            RootFind::Secant => {
                unimplemented!()
            }
        }
    }

    pub fn find_root(&mut self) -> Result<f64, RootError> {
        for i in 0 .. self.times {
            self.update();
            match self.find {
                RootBool::Find => {
                    println!("{}", i+1);
                    return Ok(self.root);
                }
                RootBool::Error => {
                    return Err(RootError::NaNRoot);
                }
                _ => (),
            }
        }
        Err(RootError::TimesUp)
    }
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
