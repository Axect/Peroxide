use anyhow::{Result, bail};

pub trait RootFindingProblem<T> {
    fn function(&self, x: f64) -> Result<f64>;
    fn initial_guess(&self) -> T;
    fn derivative(&self, x: f64) -> f64 {
        unimplemented!()
    }
    fn hessian(&self, x: f64) -> f64 {
        unimplemented!()
    }
}

pub trait RootFindingMethod<T> {
    fn step<P: RootFindingProblem<T>>(&self, problem: &P, state: T) -> Result<T>;
}

pub trait RootSolver {
    fn solve<P, F, T>(&self, problem: &P, finder: &F) -> Result<f64>
    where
        P: RootFindingProblem<T>,
        F: RootFindingMethod<T>;
}

// ┌─────────────────────────────────────────────────────────┐
//  Bisection method
// └─────────────────────────────────────────────────────────┘
pub struct BisectionMethod {
    pub max_iterations: usize,
    pub tolerance: f64,
}

impl RootFindingMethod<(f64, f64)> for BisectionMethod {
    fn step<P: RootFindingProblem<(f64, f64)>>(
        &self,
        problem: &P,
        state: (f64, f64),
    ) -> Result<(f64, f64)> {
        let (a, b) = state;
        let c = (a + b) / 2.0;

        let fa = problem.function(a)?;
        let fb = problem.function(b)?;
        let fc = problem.function(c)?;

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
pub struct NewtonMethod {
    pub max_iterations: usize,
    pub tolerance: f64,
}

impl RootFindingMethod<f64> for NewtonMethod {
    fn step<P: RootFindingProblem<f64>>(
        &self,
        problem: &P,
        state: f64,
    ) -> Result<f64> {
        let f = problem.function(state)?;
        let df = problem.derivative(state);

        if df == 0.0 {
            bail!("Zero derivative at x = {}", state);
        }

        Ok(state - f / df)
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Secant method
// └─────────────────────────────────────────────────────────┘
pub struct SecantMethod {
    pub max_iterations: usize,
    pub tolerance: f64,
}

impl RootFindingMethod<(f64, f64)> for SecantMethod {
    fn step<P: RootFindingProblem<(f64, f64)>>(
        &self,
        problem: &P,
        state: (f64, f64),
    ) -> Result<(f64, f64)> {
        let (x0, x1) = state;
        let f0 = problem.function(x0)?;
        let f1 = problem.function(x1)?;

        if f0 == f1 {
            bail!("Zero secant at ({}, {})", x0, x1);
        }

        Ok((x1, x1 - f1 * (x1 - x0) / (f1 - f0)))
    }
}
