use anyhow::{Result, bail};

type PT<const N: usize> = [f64; N];
type INTV<const N: usize> = (PT<N>, PT<N>);
type J<const R: usize, const C: usize> = [[f64; C]; R];
type H<const R: usize, const C: usize> = [[[f64; C]; C]; R];

/// Trait to define a root finding problem
///
/// # Type Parameters
///
/// - `I`: Input type (e.g. `f64`, `[f64; N]`, or etc.)
/// - `O`: Output type (e.g. `f64`, `[f64; N]`, or etc.)
/// - `T`: State type (e.g. `f64`, `(f64, f64)`, or etc.)
pub trait RootFindingProblem<const I: usize, const O: usize, T> {
    fn function(&self, x: PT<I>) -> Result<PT<O>>;
    fn initial_guess(&self) -> T;
    fn derivative(&self, x: PT<I>) -> Result<J<O, I>> {
        unimplemented!()
    }
    fn hessian(&self, x: PT<I>) -> Result<H<O, I>> {
        unimplemented!()
    }
}

pub trait RootFinder<const I: usize, const O: usize, T> {
    fn max_iter(&self) -> usize;
    fn tol(&self) -> f64;
    fn find<P: RootFindingProblem<I, O, T>>(&self, problem: &P) -> Result<PT<I>>;
}

#[derive(Debug, Copy, Clone)]
pub enum RootError<const I: usize> {
    NotConverge(PT<I>),
    NoRoot,
    ZeroDerivative(PT<I>),
}

impl std::fmt::Display for RootError<1> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RootError::NoRoot => write!(f, "There is no root in the interval"),
            RootError::NotConverge(a) => write!(f, "Not yet converge. Our guess is {:?}", a),
            RootError::ZeroDerivative(a) => write!(f, "Zero derivative in {:?}", a),
        }
    }
}

/// Macro for single function
///
/// # Description
///
/// For I=1, O=1, it is bother to write below code.
///
/// ```no_run
/// let fx = problem.function([x])?[0];
/// ```
///
/// This macro solve this problem as follows.
///
/// ```no_run
/// let fx = single_function!(problem, x);
/// ```
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
/// ```no_run
/// let fx = problem.derivative([x])?[0][0];
/// ```
///
/// This macro solve this problem as follows.
///
/// ```no_run
/// let fx = single_derivative!(problem, x);
/// ```
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

        if fa * fb > 0.0 {
            bail!(RootError::NoRoot);
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
                bail!(RootError::NoRoot);
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

#[derive(Debug, Copy, Clone)]
pub enum NewtonError {
    ZeroDerivative(f64),
    NotConverge(f64),
}

impl std::fmt::Display for NewtonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NewtonError::ZeroDerivative(x) => {
                write!(f, "Zero derivative at x = {}", x)
            }
            NewtonError::NotConverge(x) => {
                write!(f, "Not converge at x = {}", x)
            }
        }
    }
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
                bail!(NewtonError::ZeroDerivative(x));
            } else {
                x -= f / df;
            }
        }
        bail!(NewtonError::NotConverge(x));
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Secant method
// └─────────────────────────────────────────────────────────┘
/// Secant method
///
/// # Arguments
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: Absolute tolerance
///
/// # Caution
///
/// - The function should be differentiable
pub struct SecantMethod {
    pub max_iter: usize,
    pub tol: f64,
}

#[derive(Debug, Copy, Clone)]
pub enum SecantError {
    ZeroSecant(f64, f64),
    NotConverge(f64, f64),
}

impl std::fmt::Display for SecantError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SecantError::ZeroSecant(x0, x1) => {
                write!(f, "Zero secant at ({}, {})", x0, x1)
            }
            SecantError::NotConverge(x0, x1) => {
                write!(f, "Not converge at ({}, {})", x0, x1)
            }
        }
    }
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
        let mut state = problem.initial_guess();

        for _ in 0..self.max_iter {
            let (x0, x1) = state;
            let f0 = single_function!(problem, x0);
            let f1 = single_function!(problem, x1);

            if f1.abs() < self.tol {
                return Ok([x1]);
            }

            if f0 == f1 {
                bail!(SecantError::ZeroSecant(x0, x1));
            }

            state = (x1, x1 - f1 * (x1 - x0) / (f1 - f0))
        }
        let (x0, x1) = state;
        bail!(SecantError::NotConverge(x0, x1));
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  False position method
// └─────────────────────────────────────────────────────────┘
/// False position method
///
/// # Arguments
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: Absolute tolerance
///
/// # Caution
///
/// - The function should be differentiable
pub struct FalsePositionMethod {
    pub max_iter: usize,
    pub tol: f64,
}

#[derive(Debug, Copy, Clone)]
pub enum FalsePositionError {
    NoRoot,
    NotConverge(f64, f64),
}

impl std::fmt::Display for FalsePositionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FalsePositionError::NoRoot => write!(f, "There is no root in the interval"),
            FalsePositionError::NotConverge(a, b) => {
                write!(f, "Not converge in [{}, {}]", a, b)
            }
        }
    }
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

        if fa * fb > 0.0 {
            bail!(FalsePositionError::NoRoot);
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
                bail!(FalsePositionError::NoRoot);
            }
        }
        bail!(FalsePositionError::NotConverge(a, b));
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Broyden method
// └─────────────────────────────────────────────────────────┘
/// Broyden method
///
/// # Arguments
///
/// - `max_iter`: Maximum number of iterations
/// - `tol`: Absolute tolerance
///
/// # Caution
///
/// - The function should be differentiable
pub struct BroydenMethod {
    pub max_iter: usize,
    pub tol: f64,
}


impl<const I: usize, const O: usize> RootFinder<I, O, INTV<I>> for BroydenMethod {
    fn max_iter(&self) -> usize {
        self.max_iter
    }
    fn tol(&self) -> f64 {
        self.tol
    }
    fn find<P: RootFindingProblem<I, O, INTV<I>>>(
        &self,
        problem: &P,
    ) -> Result<PT<I>> {
        unimplemented!()
    }
}
