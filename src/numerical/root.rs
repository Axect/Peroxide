use thiserror::Error;

#[derive(Error, Debug, Clone)]
pub enum RootFindingError {
    #[error("No available root")]
    NoAvailableRoot,
    #[error("Zero derivative at {0}")]
    ZeroDerivative(f64),
    #[error("Constrain violated at {0}")]
    ConstraintViolation(f64),
}

pub trait RootFindingProblem<T> {
    fn function(&self, x: f64) -> Result<f64, RootFindingError>;
    fn initial_guess(&self) -> T;
    fn derivative(&self, x: f64) -> f64 {
        unimplemented!()
    }
    fn hessian(&self, x: f64) -> f64 {
        unimplemented!()
    }
}

pub trait RootFindingMethod<T> {
    fn step<P: RootFindingProblem<T>>(&self, problem: &P, state: T) -> Result<T, RootFindingError>;
}

pub trait RootSolver {
    fn solve<P, F, T>(&self, problem: &P, finder: &F) -> Result<f64, RootFindingError>
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
    ) -> Result<(f64, f64), RootFindingError> {
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
            Err(RootFindingError::NoAvailableRoot)
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
    ) -> Result<f64, RootFindingError> {
        let f = problem.function(state)?;
        let df = problem.derivative(state);

        if df == 0.0 {
            return Err(RootFindingError::ZeroDerivative(state));
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
    ) -> Result<(f64, f64), RootFindingError> {
        let (x0, x1) = state;
        let f0 = problem.function(x0)?;
        let f1 = problem.function(x1)?;

        if f0 == f1 {
            return Err(RootFindingError::ZeroDerivative(x0));
        }

        unimplemented!()
    }
}
