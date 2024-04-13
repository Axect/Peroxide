use anyhow::{Result, bail};

/// Point alias (`[f64; N]`)
pub type Pt<const N: usize> = [f64; N];
/// Interval alias (`([f64; N], [f64; N])`)
pub type Intv<const N: usize> = (Pt<N>, Pt<N>);
/// Jacobian alias (`[[f64; C]; R]`)
pub type Jaco<const R: usize, const C: usize> = [[f64; C]; R];
/// Hessian alias (`[[[f64; C]; C]; R]`)
pub type Hess<const R: usize, const C: usize> = [[[f64; C]; C]; R];

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

impl std::fmt::Display for RootError<1> {
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

            f0 = f1;
            (x0, x1) = (x1, x1 - f1 * (x1 - x0) / (f1 - f0))
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
            bail!(RootError::NoRoot);
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
                bail!(RootError::NoRoot);
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
pub struct BroydenMethod {
    pub max_iter: usize,
    pub tol: f64,
}


#[allow(unused_variables)]
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
        unimplemented!()
    }
}
