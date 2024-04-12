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

pub trait RootFindingMethod<const I: usize, const O: usize, T> {
    fn step<P: RootFindingProblem<I, O, T>>(&self, problem: &P, state: T) -> Result<T>;
}

pub trait RootSolver {
    fn solve<P, F, const I: usize, const O: usize, T>(&self, problem: &P, finder: &F) -> Result<[f64; I]>
    where
        P: RootFindingProblem<I, O, T>,
        F: RootFindingMethod<I, O, T>;
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
/// - `max_iterations`: Maximum number of iterations
/// - `tolerance`: Tolerance
///
/// # Caution
///
/// - The function should be continuous
/// - The function should have a root in the initial interval
pub struct BisectionMethod {
    pub max_iterations: usize,
    pub tolerance: f64,
}

impl RootFindingMethod<1, 1, (f64, f64)> for BisectionMethod {
    fn step<P: RootFindingProblem<1, 1, (f64, f64)>>(
        &self,
        problem: &P,
        state: (f64, f64),
    ) -> Result<(f64, f64)> {
        let (a, b) = state;
        let c = (a + b) / 2.0;

        let fa = single_function!(problem, a);
        let fb = single_function!(problem, b);
        let fc = single_function!(problem, c);

        if fa * fc < 0.0 {
            Ok((a, c))
        } else if fb * fc < 0.0 {
            Ok((c, b))
        } else if fc == 0.0 {
            Ok((c, c))
        } else {
            bail!("There is no root in the interval [{}, {}]", a, b);
        }
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Newton method
// └─────────────────────────────────────────────────────────┘
/// Newton method
///
/// # Arguments
///
/// - `max_iterations`: Maximum number of iterations
/// - `tolerance`: Tolerance
///
/// # Caution
///
/// - The function should be differentiable
/// - This method highly depends on the initial guess
/// - This method is not guaranteed to converge
pub struct NewtonMethod {
    pub max_iterations: usize,
    pub tolerance: f64,
}

impl RootFindingMethod<1, 1, f64> for NewtonMethod {
    fn step<P: RootFindingProblem<1, 1, f64>>(
        &self,
        problem: &P,
        x: f64,
    ) -> Result<f64> {
        let f = single_function!(problem, x);
        let df = single_derivative!(problem, x);

        if df == 0.0 {
            bail!("Zero derivative at x = {}", x);
        }

        Ok(x - f / df)
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Secant method
// └─────────────────────────────────────────────────────────┘
/// Secant method
///
/// # Arguments
///
/// - `max_iterations`: Maximum number of iterations
/// - `tolerance`: Tolerance
///
/// # Caution
///
/// - The function should be differentiable
pub struct SecantMethod {
    pub max_iterations: usize,
    pub tolerance: f64,
}

impl RootFindingMethod<1, 1, (f64, f64)> for SecantMethod {
    fn step<P: RootFindingProblem<1, 1, (f64, f64)>>(
        &self,
        problem: &P,
        state: (f64, f64),
    ) -> Result<(f64, f64)> {
        let (x0, x1) = state;
        let f0 = single_function!(problem, x0);
        let f1 = single_function!(problem, x1);

        if f0 == f1 {
            bail!("Zero secant at ({}, {})", x0, x1);
        }

        Ok((x1, x1 - f1 * (x1 - x0) / (f1 - f0)))
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Broyden method
// └─────────────────────────────────────────────────────────┘
/// Broyden method
///
/// # Arguments
///
/// - `max_iterations`: Maximum number of iterations
/// - `tolerance`: Tolerance
///
/// # Caution
///
/// - The function should be differentiable
pub struct BroydenMethod {
    pub max_iterations: usize,
    pub tolerance: f64,
}


impl<const I: usize, const O: usize> RootFindingMethod<I, O, INTV<I>> for BroydenMethod {
    fn step<P: RootFindingProblem<I, O, INTV<I>>>(
        &self,
        problem: &P,
        state: INTV<I>,
    ) -> Result<INTV<I>> {
        unimplemented!()
    }
}
