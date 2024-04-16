//! # Root Finding Methods
//!
//! This module provides a collection of root finding algorithms for solving nonlinear equations.
//! It defines traits for representing root finding problems and root finding methods, and provides implementations of several common algorithms.
//!
//! ## Traits
//!
//! - `RootFindingProblem<const I: usize, const O: usize, T>`: Defines the interface for a root finding problem.
//!   It requires implementing the `function` method to evaluate the function at a given point, and the `initial_guess` method to provide an initial guess for the root.
//!   Optionally, the `derivative` and `hessian` methods can be implemented to provide the derivative and Hessian of the function, respectively.
//!
//!   - `I`: Input dimension
//!   - `O`: Output dimension
//!   - `T`: State type (e.g. `f64`, `(f64, f64)`, or etc.)
//!
//! - `RootFinder<const I: usize, const O: usize, T>`: Defines the interface for a root finding method.
//!   It requires implementing the `find` method, which takes a `RootFindingProblem` and returns the root of the function.
//!   The `max_iter` and `tol` methods provide the maximum number of iterations and the tolerance for the root finding algorithm.
//!
//! ## Root Finding Methods
//!
//! - `BisectionMethod`: Implements the bisection method for finding roots of continuous functions.
//!   It requires an initial interval that brackets the root.
//!
//!  - Type Parameters: `I=1, O=1, T=(f64, f64)`
//!
//! - `NewtonMethod`: Implements Newton's method for finding roots of differentiable functions.
//!   It requires an initial guess for the root and the derivative of the function.
//!
//!   - Type Parameters: `I=1, O=1, T=f64`
//!
//! - `SecantMethod`: Implements the secant method for finding roots of differentiable functions.
//!   It requires two initial guesses for the root.
//!
//!   - Type Parameters: `I=1, O=1, T=f64`
//!
//! - `FalsePositionMethod`: Implements the false position method for finding roots of continuous functions.
//!   It requires an initial interval that brackets the root.
//!
//!   - Type Parameters: `I=1, O=1, T=(f64, f64)`
//!
//! - `BroydenMethod`: Implements Broyden's method for finding roots of systems of nonlinear equations.
//!   It requires an two initial guesses for the first step. (not an interval, just two points)
//!
//!   - Type Parameters: `I>=1, O>=1, T=([f64; I], [f64; I])`
//!
//! ## Convenient type aliases
//!
//! - `Pt<const N: usize>`: Represents a point in N-dimensional space. (`[f64; N]`)
//! - `Intv<const N: usize>`: Represents an interval in I-dimensional space. (`([f64; N], [f64; N])`)
//! - `Jaco<const R: usize, const C: usize>`: Represents the Jacobian matrix of a function. (`[[f64; C]; R]`)
//! - `Hess<const R: usize, const C: usize>`: Represents the Hessian matrix of a function. (`[[[f64; C]; C]; R]`)
//!
//! ## High-level macros
//!
//! Peroxide also provides high-level macros for root finding.
//! Assume `f: fn(f64) -> f64`.
//!
//! - `bisection!(f, (a,b), max_iter, tol)`
//! - `newton!(f, x0, max_iter, tol)`: (**Caution**: newton macro requires `#[ad_function]` attribute)
//! - `secant!(f, (x0, x1), max_iter, tol)`
//! - `false_position!(f, (a,b), max_iter, tol)`
//!
//! ```rust
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//! use anyhow::Result;
//! 
//! fn main() -> Result<()> {
//!     let root_bisect = bisection!(f, (0.0, 2.0), 100, 1e-6);
//!     let root_newton = newton!(f, 0.0, 100, 1e-6);
//!     let root_false_pos = false_position!(f, (0.0, 2.0), 100, 1e-6);
//!     let root_secant = secant!(f, (0.0, 2.0), 100, 1e-6);
//! 
//!     println!("root_bisect: {}", root_bisect);
//!     println!("root_newton: {}", root_newton);
//!     println!("root_false_pos: {}", root_false_pos);
//!     println!("root_secant: {}", root_secant);
//! 
//!     assert!(f(root_bisect).abs() < 1e-6);
//!     assert!(f(root_newton).abs() < 1e-6);
//!     assert!(f(root_false_pos).abs() < 1e-6);
//!     assert!(f(root_secant).abs() < 1e-6);
//! 
//!     Ok(())
//! }
//! 
//! #[ad_function]
//! fn f(x: f64) -> f64 {
//!     (x - 1f64).powi(3)
//! }
//! ```
//!
//! ## Examples
//!
//! ### Finding the root of a cubic function
//!
//! ```rust
//! use peroxide::fuga::*;
//! use anyhow::Result;
//!
//! fn main() -> Result<()> {
//!     let problem = Cubic;
//!
//!     let bisect = BisectionMethod { max_iter: 100, tol: 1e-6 };
//!     let newton = NewtonMethod { max_iter: 100, tol: 1e-6 };
//!     let false_pos = FalsePositionMethod { max_iter: 100, tol: 1e-6 };
//!
//!     let root_bisect = bisect.find(&problem)?;
//!     let root_newton = newton.find(&problem)?;
//!     let root_false_pos = false_pos.find(&problem)?;
//!
//!     let result_bisect = problem.eval(root_bisect)?[0];
//!     let result_newton = problem.eval(root_newton)?[0];
//!     let result_false_pos = problem.eval(root_false_pos)?[0];
//!
//!     assert!(result_bisect.abs() < 1e-6);
//!     assert!(result_newton.abs() < 1e-6);
//!     assert!(result_false_pos.abs() < 1e-6);
//!
//!     Ok(())
//! }
//!
//! struct Cubic;
//!
//! impl Cubic {
//!     fn eval(&self, x: [f64; 1]) -> Result<[f64; 1]> {
//!         Ok([(x[0] - 1f64).powi(3)])
//!     }
//! }
//!
//! impl RootFindingProblem<1, 1, (f64, f64)> for Cubic {
//!     fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
//!         self.eval(x)
//!     }
//!
//!     fn initial_guess(&self) -> (f64, f64) {
//!         (0.0, 2.0)
//!     }
//! }
//!
//! impl RootFindingProblem<1, 1, f64> for Cubic {
//!     fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
//!         self.eval(x)
//!     }
//!
//!     fn initial_guess(&self) -> f64 {
//!         0.0
//!     }
//!
//!     fn derivative(&self, x: [f64; 1]) -> Result<Jaco<1, 1>> {
//!         Ok([[3.0 * (x[0] - 1f64).powi(2)]])
//!     }
//! }
//! ```
//!
//! This example demonstrates how to find the root of a cubic function `(x - 1)^3` using various root finding methods.
//! The `Cubic` struct implements the `RootFindingProblem` trait for both `(f64, f64)` and `f64` initial guess types, allowing the use of different root finding methods.
//!
//! ### Finding the root of the cosine function (error handling example)
//!
//! ```rust
//! use peroxide::fuga::*;
//! use anyhow::Result;
//!
//! fn main() -> Result<()> {
//!     let problem = Cosine;
//!     let newton = NewtonMethod { max_iter: 100, tol: 1e-6 };
//!
//!     let root_newton = match newton.find(&problem) {
//!         Ok(x) => x,
//!         Err(e) => {
//!             println!("{:?}", e);
//!             match e.downcast::<RootError<1>>() {
//!                 Ok(RootError::ZeroDerivative(x)) => x,
//!                 Ok(e) => panic!("ok but {:?}", e),
//!                 Err(e) => panic!("err {:?}", e),
//!             }
//!         }
//!     };
//!
//!     assert_eq!(root_newton[0], 0.0);
//!
//!     Ok(())
//! }
//!
//! struct Cosine;
//!
//! impl Cosine {
//!     fn eval(&self, x: [f64; 1]) -> Result<[f64; 1]> {
//!         Ok([x[0].cos()])
//!     }
//! }
//!
//! impl RootFindingProblem<1, 1, f64> for Cosine {
//!     fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
//!         self.eval(x)
//!     }
//!
//!     fn initial_guess(&self) -> f64 {
//!         0.0 // should fail in newton (derivative is 0)
//!     }
//!
//!     fn derivative(&self, x: [f64; 1]) -> Result<Jaco<1, 1>> {
//!         Ok([[-x[0].sin()]])
//!     }
//! }
//! ```
//!
//! This example shows how to find the root of the cosine function using Newton's method.
//! The `Cosine` struct implements the `RootFindingProblem` trait for the `f64` initial guess type.
//! The initial guess is set to `0.0`, which is a point where the derivative of the cosine function is 0.
//! This leads to the `NewtonMethod` returning a `RootError::ZeroDerivative` error, which is handled in the example.
use anyhow::{Result, bail};

use crate::traits::math::{Normed, Norm, LinearOp};
use crate::traits::sugar::{ConvToMat, VecOps};
use crate::util::non_macro::zeros;

// ┌─────────────────────────────────────────────────────────┐
//  High level macro
// └─────────────────────────────────────────────────────────┘
/// High level macro for bisection
///
/// # Arguments
///
/// - `f`: `fn(f64) -> f64`
/// - `(a, b)`: `(f64, f64)`
/// - `max_iter`: `usize`
/// - `tol`: `f64`
#[macro_export]
macro_rules! bisection {
    ($f:ident, ($a:expr, $b:expr), $max_iter:expr, $tol:expr) => {{
        struct BisectionProblem;

        impl RootFindingProblem<1, 1, (f64, f64)> for BisectionProblem {
            fn initial_guess(&self) -> (f64, f64) {
                ($a, $b)
            }

            fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
                Ok([$f(x[0])])
            }
        }

        let problem = BisectionProblem;
        let bisection = BisectionMethod { max_iter: $max_iter, tol: $tol };
        let root = bisection.find(&problem)?;
        root[0]
    }}
}

/// High level macro for newton (using Automatic differentiation)
///
/// # Requirements
///
/// - This macro requires the function with `ad_function`
///
///   ```rust
///   use peroxide::fuga::*;
///
///   #[ad_function]
///   fn f(x: f64) -> f64 {
///       (x - 1f64).powi(3)
///   }
///   ```
///
/// # Arguments
///
/// - `f`: `fn(f64) -> f64`
/// - `x`: `f64`
/// - `max_iter`: `usize`
/// - `tol`: `f64`
#[macro_export]
macro_rules! newton {
    ($f:ident, $x:expr, $max_iter:expr, $tol:expr) => {{
        use paste::paste;
        struct NewtonProblem;

        impl RootFindingProblem<1, 1, f64> for NewtonProblem {
            fn initial_guess(&self) -> f64 {
                $x 
            }

            fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
                Ok([$f(x[0])])
            }

            fn derivative(&self, x: [f64; 1]) -> Result<Jaco<1, 1>> {
                paste! {
                    let x_ad = AD1(x[0], 1f64);
                    Ok([[[<$f _ad>](x_ad).dx()]])
                }
            }
        }

        let problem = NewtonProblem;
        let newton = NewtonMethod { max_iter: $max_iter, tol: $tol };
        let root = newton.find(&problem)?;
        root[0]
    }}
}

/// High level macro for false position
///
/// # Arguments
///
/// - `f`: `fn(f64) -> f64`
/// - `(a, b)`: `(f64, f64)`
/// - `max_iter`: `usize`
/// - `tol`: `f64`
#[macro_export]
macro_rules! false_position {
    ($f:ident, ($a:expr, $b:expr), $max_iter:expr, $tol:expr) => {{
        struct FalsePositionProblem;

        impl RootFindingProblem<1, 1, (f64, f64)> for FalsePositionProblem {
            fn initial_guess(&self) -> (f64, f64) {
                ($a, $b)
            }

            fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
                Ok([$f(x[0])])
            }
        }

        let problem = FalsePositionProblem;
        let false_position = FalsePositionMethod { max_iter: $max_iter, tol: $tol };
        let root = false_position.find(&problem)?;
        root[0]
    }}
}

/// High level macro for secant
///
/// # Arguments
///
/// - `f`: `fn(f64) -> f64`
/// - `(a, b)`: `(f64, f64)`
/// - `max_iter`: `usize`
/// - `tol`: `f64`
#[macro_export]
macro_rules! secant {
    ($f:ident, ($a:expr, $b:expr), $max_iter:expr, $tol:expr) => {{
        struct SecantProblem;

        impl RootFindingProblem<1, 1, (f64, f64)> for SecantProblem {
            fn initial_guess(&self) -> (f64, f64) {
                ($a, $b)
            }

            fn function(&self, x: [f64; 1]) -> Result<[f64; 1]> {
                Ok([$f(x[0])])
            }
        }

        let problem = SecantProblem;
        let secant = SecantMethod { max_iter: $max_iter, tol: $tol };
        let root = secant.find(&problem)?;
        root[0]
    }}
}


// ┌─────────────────────────────────────────────────────────┐
//  Type aliases
// └─────────────────────────────────────────────────────────┘
/// Point alias (`[f64; N]`)
pub type Pt<const N: usize> = [f64; N];
/// Interval alias (`([f64; N], [f64; N])`)
pub type Intv<const N: usize> = (Pt<N>, Pt<N>);
/// Jacobian alias (`[[f64; C]; R]`)
pub type Jaco<const R: usize, const C: usize> = [[f64; C]; R];
/// Hessian alias (`[[[f64; C]; C]; R]`)
pub type Hess<const R: usize, const C: usize> = [[[f64; C]; C]; R];

// ┌─────────────────────────────────────────────────────────┐
//  Traits
// └─────────────────────────────────────────────────────────┘
/// Trait to define a root finding problem
///
/// # Type Parameters
///
/// - `I`: Input type (e.g. `f64`, `[f64; N]`, or etc.)
/// - `O`: Output type (e.g. `f64`, `[f64; N]`, or etc.)
/// - `T`: State type (e.g. `f64`, `(f64, f64)`, or etc.)
///
/// # Methods
///
/// - `function`: Function
/// - `initial_guess`: Initial guess
/// - `derivative`: Derivative (optional)
/// - `hessian`: Hessian (optional)
pub trait RootFindingProblem<const I: usize, const O: usize, T> {
    fn function(&self, x: Pt<I>) -> Result<Pt<O>>;
    fn initial_guess(&self) -> T;
    #[allow(unused_variables)]
    fn derivative(&self, x: Pt<I>) -> Result<Jaco<O, I>> {
        unimplemented!()
    }
    #[allow(unused_variables)]
    fn hessian(&self, x: Pt<I>) -> Result<Hess<O, I>> {
        unimplemented!()
    }
}

/// Trait to define a root finder
///
/// # Type Parameters
///
/// - `I`: Input type (e.g. `f64`, `[f64; N]`, or etc.)
/// - `O`: Output type (e.g. `f64`, `[f64; N]`, or etc.)
/// - `T`: State type (e.g. `f64`, `(f64, f64)`, or etc.)
///
/// # Methods
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: Absolute tolerance
/// - `find`: Find root
///
/// # Available root finders
///
/// - `BisectionMethod`: `I=1, O=1, T=(f64, f64)`
/// - `FalsePositionMethod`: `I=1, O=1, T=(f64, f64)`
/// - `NewtonMethod`: `I=1, O=1, T=f64`
/// - `SecantMethod`: `I=1, O=1, T=(f64, f64)`
pub trait RootFinder<const I: usize, const O: usize, T> {
    fn max_iter(&self) -> usize;
    fn tol(&self) -> f64;
    fn find<P: RootFindingProblem<I, O, T>>(&self, problem: &P) -> Result<Pt<I>>;
}

#[derive(Debug, Copy, Clone)]
pub enum RootError<const I: usize> {
    NotConverge(Pt<I>),
    NoRoot,
    ZeroDerivative(Pt<I>),
    ZeroSecant(Pt<I>, Pt<I>),
}

impl<const I: usize> std::fmt::Display for RootError<I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RootError::NoRoot => write!(f, "There is no root in the interval"),
            RootError::NotConverge(a) => write!(f, "Not yet converge. Our guess is {:?}", a),
            RootError::ZeroDerivative(a) => write!(f, "Zero derivative in {:?}", a),
            RootError::ZeroSecant(a,b) => write!(f, "Zero secant in ({:?}, {:?})", a, b),
        }
    }
}

/// Macro for single function
///
/// # Description
///
/// For I=1, O=1, it is bother to write below code.
///
/// ```ignore
/// let fx = problem.function([x])?[0];
/// ```
///
/// This macro solve this problem as follows.
///
/// ```ignore
/// let fx = single_function!(problem, x);
/// ```
#[macro_export]
macro_rules! single_function {
    ($problem:expr, $x:expr) => {{
        $problem.function([$x])?[0]
    }}
}

/// Macro for single derivative
///
/// # Description
///
/// For I=1, O=1, it is bother to write below code.
///
/// ```ignore
/// let fx = problem.derivative([x])?[0][0];
/// ```
///
/// This macro solve this problem as follows.
///
/// ```ignore
/// let fx = single_derivative!(problem, x);
/// ```
#[macro_export]
macro_rules! single_derivative {
    ($problem:expr, $x:expr) => {{
        $problem.derivative([$x])?[0][0]
    }}
}

// ┌─────────────────────────────────────────────────────────┐
//  Bisection method
// └─────────────────────────────────────────────────────────┘
/// Bisection method
///
/// # Type for `RootFinder`
///
/// - `I`: 1
/// - `O`: 1
/// - `T`: `(f64, f64)`
///
/// # Arguments
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: tol
///
/// # Caution
///
/// - The function should be continuous
/// - The function should have a root in the initial interval
pub struct BisectionMethod {
    pub max_iter: usize,
    pub tol: f64,
}

impl RootFinder<1, 1, (f64, f64)> for BisectionMethod {
    fn max_iter(&self) -> usize {
        self.max_iter
    }

    fn tol(&self) -> f64 {
        self.tol
    }

    fn find<P: RootFindingProblem<1, 1, (f64, f64)>>(
        &self,
        problem: &P,
    ) -> Result<[f64; 1]> {
        let state = problem.initial_guess();
        let (mut a, mut b) = state;
        let mut fa = single_function!(problem, a);
        let mut fb = single_function!(problem, b);

        if fa.abs() < self.tol {
            return Ok([a]);
        } else if fb.abs() < self.tol {
            return Ok([b]);
        } else if fa * fb > 0.0 {
            bail!(RootError::<1>::NoRoot);
        }

        for _ in 0..self.max_iter {
            let c = (a + b) / 2.0;
            let fc = single_function!(problem, c);

            if fc.abs() < self.tol {
                return Ok([c]);
            } else if fa * fc < 0.0 {
                b = c;
                fb = fc;
            } else if fb * fc < 0.0 {
                a = c;
                fa = fc;
            } else {
                bail!(RootError::<1>::NoRoot);
            }
        }
        let c = (a + b) / 2.0;
        bail!(RootError::NotConverge([c]));
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Newton method
// └─────────────────────────────────────────────────────────┘
/// Newton method
///
/// # Type for `RootFinder`
///
/// - `I`: 1
/// - `O`: 1
/// - `T`: `f64`
///
/// # Arguments
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: Absolute tolerance
///
/// # Caution
///
/// - The function should be differentiable
/// - This method highly depends on the initial guess
/// - This method is not guaranteed to converge
pub struct NewtonMethod {
    pub max_iter: usize,
    pub tol: f64,
}

impl RootFinder<1, 1, f64> for NewtonMethod {
    fn max_iter(&self) -> usize {
        self.max_iter
    }
    fn tol(&self) -> f64 {
        self.tol
    }
    fn find<P: RootFindingProblem<1, 1, f64>>(
        &self,
        problem: &P,
    ) -> Result<[f64; 1]> {
        let mut x = problem.initial_guess();

        for _ in 0..self.max_iter {
            let f = single_function!(problem, x);
            let df = single_derivative!(problem, x);

            if f.abs() < self.tol {
                return Ok([x]);
            } else if df == 0.0 {
                bail!(RootError::ZeroDerivative([x]));
            } else {
                x -= f / df;
            }
        }
        bail!(RootError::NotConverge([x]));
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Secant method
// └─────────────────────────────────────────────────────────┘
/// Secant method
///
/// # Type for `RootFinder`
///
/// - `I`: 1
/// - `O`: 1
/// - `T`: `(f64, f64)`
///
/// # Arguments
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: Absolute tolerance
///
/// # Caution
///
/// - The function should be differentiable
/// - This method is not guaranteed to converge
pub struct SecantMethod {
    pub max_iter: usize,
    pub tol: f64,
}

impl RootFinder<1, 1, (f64, f64)> for SecantMethod {
    fn max_iter(&self) -> usize {
        self.max_iter
    }
    fn tol(&self) -> f64 {
        self.tol
    }
    fn find<P: RootFindingProblem<1, 1, (f64, f64)>>(
        &self,
        problem: &P,
    ) -> Result<[f64; 1]> {
        let state = problem.initial_guess();
        let (mut x0, mut x1) = state;
        let mut f0 = single_function!(problem, x0);

        if f0.abs() < self.tol {
            return Ok([x0]);
        }

        for _ in 0..self.max_iter {
            let f1 = single_function!(problem, x1);

            if f1.abs() < self.tol {
                return Ok([x1]);
            }

            if f0 == f1 {
                bail!(RootError::ZeroSecant([x0], [x1]));
            }

            let f0_old = f0;
            f0 = f1;
            (x0, x1) = (x1, x1 - f1 * (x1 - x0) / (f1 - f0_old))
        }
        bail!(RootError::NotConverge([x1]));
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  False position method
// └─────────────────────────────────────────────────────────┘
/// False position method
///
/// # Type for `RootFinder`
///
/// - `I`: 1
/// - `O`: 1
/// - `T`: `(f64, f64)`
///
/// # Arguments
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: Absolute tolerance
///
/// # Caution
///
/// - The function should be continuous
pub struct FalsePositionMethod {
    pub max_iter: usize,
    pub tol: f64,
}

impl RootFinder<1, 1, (f64, f64)> for FalsePositionMethod {
    fn max_iter(&self) -> usize {
        self.max_iter
    }
    fn tol(&self) -> f64 {
        self.tol
    }
    fn find<P: RootFindingProblem<1, 1, (f64, f64)>>(
        &self,
        problem: &P,
    ) -> Result<[f64; 1]> {
        let state = problem.initial_guess();
        let (mut a, mut b) = state;
        let mut fa = single_function!(problem, a);
        let mut fb = single_function!(problem, b);

        if fa.abs() < self.tol {
            return Ok([a]);
        } else if fb.abs() < self.tol {
            return Ok([b]);
        } else if fa * fb > 0.0 {
            bail!(RootError::<1>::NoRoot);
        }

        for _ in 0..self.max_iter {
            let c = (a * fb - b * fa) / (fb - fa);
            let fc = single_function!(problem, c);

            if fc.abs() < self.tol {
                return Ok([c]);
            } else if fa * fc < 0.0 {
                b = c;
                fb = fc;
            } else if fb * fc < 0.0 {
                a = c;
                fa = fc;
            } else {
                bail!(RootError::<1>::NoRoot);
            }
        }
        let c = (a * fb - b * fa) / (fb - fa);
        bail!(RootError::NotConverge([c]));
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Broyden method
// └─────────────────────────────────────────────────────────┘
/// Broyden method
///
/// # Type for `RootFinder`
///
/// - `I`: free
/// - `O`: free
/// - `T`: `Intv<I>` (=`([f64; I], [f64; I])`)
///
/// # Arguments
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: Absolute tolerance
///
/// # Caution
///
/// - The function should be differentiable
///
/// # Example
///
/// ```rust
/// use peroxide::fuga::*;
/// use peroxide::numerical::root::{Pt, Intv};
/// 
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let problem = CircleTangentLine;
///     let broyden = BroydenMethod { max_iter: 100, tol: 1e-6 };
///
///     let root = broyden.find(&problem)?;
///     let result = problem.function(root)?;
///
///     let norm = result.to_vec().norm(Norm::L2);
///     assert!(norm < 1e-6);
///
///     Ok(())
/// }
/// 
/// struct CircleTangentLine;
/// 
/// impl RootFindingProblem<2, 2, Intv<2>> for CircleTangentLine {
///     fn function(&self, x: Pt<2>) -> anyhow::Result<Pt<2>> {
///         Ok([
///             x[0] * x[0] + x[1] * x[1] - 1.0,
///             x[0] + x[1] - 2f64.sqrt()
///         ])
///     }
/// 
///     fn initial_guess(&self) -> Intv<2> {
///         ([0.0, 0.1], [-0.1, 0.2])
///     }
/// }
/// ```

pub struct BroydenMethod {
    pub max_iter: usize,
    pub tol: f64,
    pub rtol: f64,
}

#[allow(unused_variables, non_snake_case)]
impl<const I: usize, const O: usize> RootFinder<I, O, Intv<I>> for BroydenMethod {
    fn max_iter(&self) -> usize {
        self.max_iter
    }
    fn tol(&self) -> f64 {
        self.tol
    }
    fn find<P: RootFindingProblem<I, O, Intv<I>>>(
        &self,
        problem: &P,
    ) -> Result<Pt<I>> {
        // Init state
        let state = problem.initial_guess();
        let (mut x0, mut x1) = state;
        let mut fx0 = problem.function(x0)?.to_vec();

        // Initialize negative inverse jacobian as identity
        // H = -J^{-1}
        let mut H = zeros(I, O);
        for i in 0..O.min(I) {
            H[(i, i)] = 1.0;
        }

        for _ in 0..self.max_iter {
            let fx1 = problem.function(x1)?.to_vec();
            if fx1.norm(Norm::L2) < self.tol {
                return Ok(x1);
            }
            let dx = x1.iter().zip(x0.iter()).map(|(x1, x0)| x1 - x0).collect::<Vec<_>>();
            if dx.norm(Norm::L2) < self.rtol {
                return Ok(x1);
            }
            let df = fx1.iter().zip(fx0.iter()).map(|(fx1, fx0)| fx1 - fx0).collect::<Vec<_>>();

            let denom = dx.add_v(&H.apply(&df));
            let right = &dx.to_row() * &H;
            let num = right.apply(&df)[0];

            let left = denom.div_s(num);

            H = H - left.to_col() * right;

            x0 = x1;
            fx0 = fx1.clone();
            let dx_new = H.apply(&fx1);
            for i in 0..I {
                x1[i] += dx_new[i];
            }
        }
        bail!(RootError::NotConverge(x1));
    }
}
