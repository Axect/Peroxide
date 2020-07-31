//! Root Finding
//!
//! # Implemented Algorithms
//!
//! * Bisection
//! * False Position (Regula Falsi)
//! * Secant
//! * Newton
//!
//! # Low-level API
//!
//! ## Members
//!
//! * `RootState`
//!     * `P(f64)` : For point-like initial guess
//!     * `I(f64, f64)` : For interval-like initial guess
//! * `RootFind` : Algorithms for root finding
//!     * `Bisection`
//!     * `FalsePosition`
//!     * `Newton`
//!     * `Secant`
//! * `RootError` : Error for root finding
//!     * `MismatchedState` : Mismatched state and method (ex: Point vs Bisection)
//!     * `TimesUp` : No root until `self.times`
//!     * `NaNRoot` : NaN
//! * `RootFinder` : Main structure for root finding
//!     * `fn new(RootState, RootFind, f)` : Creat RootFinder (times: 100, tol: 1e-10)
//!     * `fn condition_number(&self) -> f64` : Compute condition number
//!     * `fn set_tol(&mut self, f64) -> &mut Self` : Set tolerance
//!     * `fn set_times(&mut self, usize) -> &mut Self` : Set max iteration times
//!     * `fn update(&mut self)` : Update one step
//!     * `fn find_root(&mut self) -> Result<f64, RootError>` : Find root
//!
//! ## Usage
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() -> Result<(), RootError> {
//!     let init = RootState::I(1f64, 4f64);
//!     let mut rf = RootFinder::<_, AD1>::new(init, Bisection, f)?;
//!     rf.set_tol(1e-15)       // Default: 1e-10
//!         .set_times(200);    // Default: 100
//!     let root = rf.find_root()?;
//!     root.print();
//!     Ok(())
//! }
//!
//! fn f<T: AD>(x: T) -> T {
//!     x.sin()
//! }
//! ```
//!
//! # High-level API
//!
//! ## Members
//!
//! All output type is `Result<f64, RootError>`
//!
//! * `bisection(f, interval: (f64, f64), times: usize, tol: f64)`
//! * `false_position(f, interval: (f64, f64), times: usize, tol: f64)`
//! * `secant(f, initial_guess: (f64, f64), times: usize, tol: f64)`
//! * `newton(f, initial_guess: f64, times: usize, tol: f64)`
//!
//! ## Usage
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() -> Result<(), RootError> {
//!     let root = bisection(f, (1f64, 4f64), 100, 1e-15)?;
//!     root.print();
//!     Ok(())
//! }
//!
//! fn f<T: AD>(x: T) -> T {
//!     x.sin()
//! }
//! ```
//!
//! # Reference
//!
//! * Walter Gautschi, *Numerical Analysis*, Springer (2012)

use crate::traits::stable::StableFn;
use crate::structure::ad::{AD, AD1, ADLift};
use std::marker::PhantomData;
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
pub struct RootFinder<F: Fn(T) -> T, T: AD> {
    init: RootState,
    pub curr: RootState, 
    method: RootFind,
    f: Box<F>,
    find: RootBool,
    times: usize,
    tol: f64,
    root: f64,
    _marker: PhantomData<T>
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

impl<F: Fn(T) -> T, T: AD> RootFinder<F, T> {
    /// Create RootFinder
    ///
    /// # Default Options
    ///
    /// * Times: 100
    /// * Tolerance: 1e-10
    ///
    /// # Usage
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() -> Result<(), RootError> {
    ///     let init = RootState::I(1f64, 3f64);
    ///     let mut a = RootFinder::<_, AD1>::new(init, Bisection, f)?;
    ///     let x = a.find_root()?;
    ///     x.print();
    ///     Ok(())
    /// }
    ///
    /// fn f<T: AD>(x: T) -> T {
    ///     x.powi(2) - 4f64
    /// }
    /// ```
    pub fn new(init: RootState, method: RootFind, f: F) -> Result<Self, RootError> {
        match method {
            RootFind::Bisection => {
                match init {
                    P(_) => Err(RootError::MismatchedState),
                    _ => Ok(
                        RootFinder {
                            init,
                            curr: init,
                            method,
                            f: Box::new(f),
                            find: RootBool::NotYet,
                            times: 100,
                            tol: 1e-10,
                            root: 0f64,
                            _marker: PhantomData
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
                            f: Box::new(f),
                            find: RootBool::NotYet,
                            times: 100,
                            tol: 1e-10,
                            root: 0f64,
                            _marker: PhantomData
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
                            f: Box::new(f),
                            find: RootBool::NotYet,
                            times: 100,
                            tol: 1e-10,
                            root: 0f64,
                            _marker: PhantomData
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
                            f: Box::new(f),
                            find: RootBool::NotYet,
                            times: 100,
                            tol: 1e-10,
                            root: 0f64,
                            _marker: PhantomData
                        }
                    ),
                    _ => Err(RootError::MismatchedState),
                }
            }
        }
    }

    pub fn f(&self, x: T) -> T {
        (self.f)(x)
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

    /// Set max iteration times
    pub fn set_times(&mut self, times: usize) -> &mut Self {
        self.times = times;
        self
    }
    
    /// Set tolerance
    pub fn set_tol(&mut self, tol: f64) -> &mut Self {
        self.tol = tol;
        self
    }

    /// Update RootFinder
    #[inline]
    pub fn update(&mut self) {
        match self.method {
            RootFind::Bisection => {
                match self.curr {
                    I(a, b) => {
                        let lift = ADLift::new(|x| self.f(x));
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
                                self.find = RootBool::Error;
                            }
                        }
                    }
                    _ => unreachable!()
                }    
            },
            RootFind::FalsePosition => {
                match self.curr {
                    I(a, b) => {
                        let lift = ADLift::new(|x| self.f(x));
                        let fa = lift.call_stable(a);
                        let fb = lift.call_stable(b);
                        let x = (a * fb - b * fa) / (fb - fa);
                        let fx = lift.call_stable(x);
                        if (a - b).abs() <= self.tol || fx.abs() <= self.tol {
                            self.find = RootBool::Find;
                            self.root = x;
                        } else {
                            if fx * fa < 0f64 {
                                self.curr = I(a, x);
                            } else if fx * fb < 0f64  {
                                self.curr = I(x, b);
                            } else if fx == 0f64 {
                                self.find = RootBool::Find;
                                self.root = x;
                            } else {
                                self.find = RootBool::Error;
                            }
                        }
                    }
                    _ => unreachable!()
                }
            },
            RootFind::Newton => {
                match self.curr {
                    P(xn) => {
                        let lift = ADLift::new(|x| self.f(x));
                        let mut z = AD1::from(xn);
                        z.d1 = 1f64;
                        let fz = lift.call_stable(z);
                        let x = xn - (fz.d0 / fz.d1);
                        if (x - xn).abs() <= self.tol {
                            self.find = RootBool::Find;
                            self.root = x;
                        }
                        self.curr = P(x);
                    },
                    _ => unreachable!()
                }
            },
            RootFind::Secant => {
                match self.curr {
                    I(xn_1, xn) => {
                        let lift = ADLift::new(|x| self.f(x));
                        let fxn_1 = lift.call_stable(xn_1);
                        let fxn = lift.call_stable(xn);
                        let x = xn - (xn - xn_1) / (fxn - fxn_1) * fxn;
                        if (x - xn).abs() <= self.tol {
                            self.find = RootBool::Find;
                            self.root = x;
                        }
                        self.curr = I(xn, x);
                    }
                    _ => unreachable!()
                }
            }
        }
    }

    /// Find Root
    pub fn find_root(&mut self) -> Result<f64, RootError> {
        for _i in 0 .. self.times {
            self.update();
            match self.find {
                RootBool::Find => {
                    //println!("{}", i+1);
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

/// Bisection method to find root
///
/// # Usage
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() -> Result<(), RootError> {
///     let x = bisection(f, (0f64, 4f64), 100, 1e-15)?;
///     assert!((x - 3f64).abs() < 1e-15);
///     Ok(())
/// }
///
/// fn f<T: AD>(x: T) -> T {
///     x.powi(2) - x * 2f64 - 3f64
/// }
/// ```
pub fn bisection<F: Fn(AD1) -> AD1>(f: F, interval: (f64, f64), times: usize, tol: f64) -> Result<f64, RootError> {
    let (a, b) = interval;
    let mut root_finder = RootFinder::new(
        RootState::I(a, b), 
        RootFind::Bisection, 
        f
    )?;
    root_finder.set_tol(tol);
    root_finder.set_times(times);
    root_finder.find_root()
}

/// False position method to find root
///
/// # Usage
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() -> Result<(), RootError> {
///     let x = false_position(f, (0f64, 4f64), 1000, 1e-15)?;
///     assert!((x - 3f64).abs() < 1e-15);
///     Ok(())
/// }
///
/// fn f<T: AD>(x: T) -> T {
///     x.powi(2) - x * 2f64 - 3f64
/// }
/// ```
pub fn false_position<F: Fn(AD1) -> AD1>(f: F, interval: (f64, f64), times: usize, tol: f64) -> Result<f64, RootError> {
    let (a, b) = interval;
    let mut root_finder = RootFinder::new(
        RootState::I(a, b), 
        RootFind::FalsePosition, 
        f
    )?;
    root_finder.set_tol(tol);
    root_finder.set_times(times);
    root_finder.find_root()
}

/// Newton method to find root
///
/// # Usage
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() -> Result<(), RootError> {
///     let x = newton(f, 2f64, 100, 1e-15)?;
///     assert!((x - 3f64).abs() < 1e-15);
///     Ok(())
/// }
///
/// fn f<T: AD>(x: T) -> T {
///     x.powi(2) - x * 2f64 - 3f64
/// }
/// ```
pub fn newton<F: Fn(AD1) -> AD1>(f: F, initial_guess: f64, times: usize, tol: f64) -> Result<f64, RootError> {
    let mut root_finder = RootFinder::new(
        RootState::P(initial_guess),
        RootFind::Newton, 
        f
    )?;
    root_finder.set_tol(tol);
    root_finder.set_times(times);
    root_finder.find_root()
}

/// Secant method to find root
///
/// # Usage
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() -> Result<(), RootError> {
///     let x = secant(f, (2f64, 2.1f64), 100, 1e-15)?;
///     assert!((x - 3f64).abs() < 1e-15);
///     Ok(())
/// }
///
/// fn f<T: AD>(x: T) -> T {
///     x.powi(2) - x * 2f64 - 3f64
/// }
/// ```
pub fn secant<F: Fn(AD1) -> AD1>(f: F, initial_guess: (f64, f64), times: usize, tol: f64) -> Result<f64, RootError> {
    let (a, b) = initial_guess;
    let mut root_finder = RootFinder::new(
        RootState::I(a, b), 
        RootFind::Secant, 
        f
    )?;
    root_finder.set_tol(tol);
    root_finder.set_times(times);
    root_finder.find_root()
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
