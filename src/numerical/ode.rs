//! # Ordinary Differential Equation (ODE) Solvers
//!
//! This module provides traits and structs for solving ordinary differential equations (ODEs).
//!
//! ## Overview
//!
//! - `ODEProblem`: Trait for defining an ODE problem.
//! - `ODEIntegrator`: Trait for ODE integrators.
//! - `ODESolver`: Trait for ODE solvers.
//! - `ODEError`: Enum for ODE errors.
//!   - `ReachedMaxStepIter`: Reached maximum number of steps per step. (internal error)
//!   - `ConstraintViolation(f64, Vec<f64>, Vec<f64>)`: Constraint violation. (user-defined error)
//!   - ODE uses `anyhow` for error handling. So, you can customize your errors.
//!
//! ## Available integrators
//!
//! - **Explicit**
//!   - Ralston's 3rd order (RALS3)
//!   - Runge-Kutta 4th order (RK4)
//!   - Ralston's 4th order (RALS4)
//!   - Runge-Kutta 5th order (RK5)
//! - **Embedded**
//!   - Bogacki-Shampine 2/3rd order (BS23)
//!   - Runge-Kutta-Fehlberg 4/5th order (RKF45)
//!   - Dormand-Prince 4/5th order (DP45)
//!   - Tsitouras 4/5th order (TSIT45)
//!   - Runge-Kutta-Fehlberg 7/8th order (RKF78)
//! - **Implicit**
//!   - Gauss-Legendre 4th order (GL4)
//!
//! ## Available solvers
//!
//! - `BasicODESolver`: A basic ODE solver using a specified integrator.
//!
//! You can implement your own ODE solver by implementing the `ODESolver` trait.
//!
//! ## Example
//!
//! ```rust
//! use peroxide::fuga::*;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // Same as : let rkf = RKF45::new(1e-4, 0.9, 1e-6, 1e-1, 100);
//!     let rkf = RKF45 {
//!         tol: 1e-6,
//!         safety_factor: 0.9,
//!         min_step_size: 1e-6,
//!         max_step_size: 1e-1,
//!         max_step_iter: 100,
//!     };
//!     let basic_ode_solver = BasicODESolver::new(rkf);
//!     let initial_conditions = vec![1f64];
//!     let (t_vec, y_vec) = basic_ode_solver.solve(
//!         &Test,
//!         (0f64, 10f64),
//!         0.01,
//!         &initial_conditions,
//!     )?;
//!     let y_vec: Vec<f64> = y_vec.into_iter().flatten().collect();
//!     println!("{}", y_vec.len());
//!
//! #   #[cfg(feature = "plot")]
//! #   {
//!     let mut plt = Plot2D::new();
//!     plt
//!         .set_domain(t_vec)
//!         .insert_image(y_vec)
//!         .set_xlabel(r"$t$")
//!         .set_ylabel(r"$y$")
//!         .set_style(PlotStyle::Nature)
//!         .tight_layout()
//!         .set_dpi(600)
//!         .set_path("example_data/rkf45_test.png")
//!         .savefig()?;
//! #   }
//!     Ok(())
//! }
//!
//! // Extremely customizable struct
//! struct Test;
//!
//! impl ODEProblem for Test {
//!     fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> anyhow::Result<()> {
//!         Ok(dy[0] = (5f64 * t.powi(2) - y[0]) / (t + y[0]).exp())
//!     }
//! }
//! ```

use crate::fuga::ConvToMat;
use crate::traits::math::{InnerProduct, Norm, Normed, Vector};
use crate::util::non_macro::eye;
use anyhow::{bail, Result};

/// Trait for defining an ODE problem.
///
/// Implement this trait to define your own ODE problem.
///
/// # Example
///
/// ```
/// use peroxide::fuga::*;
///
/// struct MyODEProblem;
///
/// impl ODEProblem for MyODEProblem {
///     fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> anyhow::Result<()> {
///         dy[0] = -0.5 * y[0];
///         dy[1] = y[0] - y[1];
///         Ok(())
///     }
/// }
/// ```
pub trait ODEProblem {
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> Result<()>;
}

/// Trait for ODE integrators.
///
/// Implement this trait to define your own ODE integrator.
pub trait ODEIntegrator {
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64>;
}

/// Enum for ODE errors.
///
/// # Variants
///
/// - `ReachedMaxStepIter`: Reached maximum number of steps per step. (internal error for integrator)
/// - `ConstraintViolation`: Constraint violation. (user-defined error)
///
/// If you define constraints in your problem, you can use this error to report constraint violations.
///
/// # Example
///
/// ```no_run
/// use peroxide::fuga::*;
///
/// struct ConstrainedProblem {
///     y_constraint: f64
/// }
///
/// impl ODEProblem for ConstrainedProblem {
///     fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> anyhow::Result<()> {
///         if y[0] < self.y_constraint {
///             anyhow::bail!(ODEError::ConstraintViolation(t, y.to_vec(), dy.to_vec()));
///         } else {
///             // some function
///             Ok(())
///         }
///     }
/// }
/// ```
#[derive(Debug, Clone)]
pub enum ODEError {
    ConstraintViolation(f64, Vec<f64>, Vec<f64>), // t, y, dy
    ReachedMaxStepIter,
}

impl std::fmt::Display for ODEError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ODEError::ConstraintViolation(t, y, dy) => write!(
                f,
                "Constraint violation at t = {}, y = {:?}, dy = {:?}",
                t, y, dy
            ),
            ODEError::ReachedMaxStepIter => write!(f, "Reached maximum number of steps per step"),
        }
    }
}

/// Trait for ODE solvers.
///
/// Implement this trait to define your own ODE solver.
pub trait ODESolver {
    fn solve<P: ODEProblem>(
        &self,
        problem: &P,
        t_span: (f64, f64),
        dt: f64,
        initial_conditions: &[f64],
    ) -> Result<(Vec<f64>, Vec<Vec<f64>>)>;
}

/// A basic ODE solver using a specified integrator.
///
/// # Example
///
/// ```
/// use peroxide::fuga::*;
///
/// fn main() -> Result<(), Box<dyn Error>> {
///     let initial_conditions = vec![1f64];
///     let rkf = RKF45::new(1e-4, 0.9, 1e-6, 1e-1, 100);
///     let basic_ode_solver = BasicODESolver::new(rkf);
///     let (t_vec, y_vec) = basic_ode_solver.solve(
///         &Test,
///         (0f64, 10f64),
///         0.01,
///         &initial_conditions,
///     )?;
///     let y_vec: Vec<f64> = y_vec.into_iter().flatten().collect();
///
///     Ok(())
/// }
///
/// struct Test;
///
/// impl ODEProblem for Test {
///     fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> anyhow::Result<()> {
///         dy[0] = (5f64 * t.powi(2) - y[0]) / (t + y[0]).exp();
///         Ok(())
///     }
/// }
/// ```
pub struct BasicODESolver<I: ODEIntegrator> {
    integrator: I,
}

impl<I: ODEIntegrator> BasicODESolver<I> {
    pub fn new(integrator: I) -> Self {
        Self { integrator }
    }
}

impl<I: ODEIntegrator> ODESolver for BasicODESolver<I> {
    fn solve<P: ODEProblem>(
        &self,
        problem: &P,
        t_span: (f64, f64),
        dt: f64,
        initial_conditions: &[f64],
    ) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
        let mut t = t_span.0;
        let mut dt = dt;
        let mut y = initial_conditions.to_vec();
        let mut t_vec = vec![t];
        let mut y_vec = vec![y.clone()];

        while t < t_span.1 {
            let dt_step = self.integrator.step(problem, t, &mut y, dt)?;
            t += dt;
            t_vec.push(t);
            y_vec.push(y.clone());
            dt = dt_step;
        }

        Ok((t_vec, y_vec))
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Butcher Tableau
// └─────────────────────────────────────────────────────────┘
/// Trait for Butcher tableau
///
/// ```text
/// C | A
/// - - -
///   | BU (Coefficient for update)
///   | BE (Coefficient for estimate error)
/// ```
///
/// # References
///
/// - J. R. Dormand and P. J. Prince, _A family of embedded Runge-Kutta formulae_, J. Comp. Appl. Math., 6(1), 19-26, 1980.
/// - Wikipedia: [List of Runge-Kutta methods](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)
pub trait ButcherTableau {
    const C: &'static [f64];
    const A: &'static [&'static [f64]];
    const BU: &'static [f64];
    const BE: &'static [f64];

    fn tol(&self) -> f64 {
        unimplemented!()
    }

    fn safety_factor(&self) -> f64 {
        unimplemented!()
    }

    fn max_step_size(&self) -> f64 {
        unimplemented!()
    }

    fn min_step_size(&self) -> f64 {
        unimplemented!()
    }

    fn max_step_iter(&self) -> usize {
        unimplemented!()
    }
}

impl<BU: ButcherTableau> ODEIntegrator for BU {
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64> {
        let n = y.len();
        let mut iter_count = 0usize;
        let mut dt = dt;
        let n_k = Self::C.len();

        loop {
            let mut k_vec = vec![vec![0.0; n]; n_k];
            let mut y_temp = y.to_vec();

            for stage in 0..n_k {
                for i in 0..n {
                    let mut s = 0.0;
                    for j in 0..stage {
                        s += Self::A[stage][j] * k_vec[j][i];
                    }
                    y_temp[i] = y[i] + dt * s;
                }
                problem.rhs(t + dt * Self::C[stage], &y_temp, &mut k_vec[stage])?;
            }

            if !Self::BE.is_empty() {
                let mut error = 0f64;
                for i in 0..n {
                    let mut s = 0.0;
                    for j in 0..n_k {
                        s += (Self::BU[j] - Self::BE[j]) * k_vec[j][i];
                    }
                    error = error.max(dt * s.abs())
                }

                let factor = (self.tol() * dt / error).powf(0.2);
                let new_dt = self.safety_factor() * dt * factor;
                let new_dt = new_dt.clamp(self.min_step_size(), self.max_step_size());

                if error < self.tol() {
                    for i in 0..n {
                        let mut s = 0.0;
                        for j in 0..n_k {
                            s += Self::BU[j] * k_vec[j][i];
                        }
                        y[i] += dt * s;
                    }
                    return Ok(new_dt);
                } else {
                    iter_count += 1;
                    if iter_count >= self.max_step_iter() {
                        bail!(ODEError::ReachedMaxStepIter);
                    }
                    dt = new_dt;
                }
            } else {
                for i in 0..n {
                    let mut s = 0.0;
                    for j in 0..n_k {
                        s += Self::BU[j] * k_vec[j][i];
                    }
                    y[i] += dt * s;
                }
                return Ok(dt);
            }
        }
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Runge-Kutta
// └─────────────────────────────────────────────────────────┘
/// Ralston's 3rd order integrator
///
/// This integrator uses the Ralston's 3rd order method to numerically integrate the ODE system.
/// In MATLAB, it is called `ode3`.
#[derive(Debug, Clone, Copy, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct RALS3;

impl ButcherTableau for RALS3 {
    const C: &'static [f64] = &[0.0, 0.5, 0.75];
    const A: &'static [&'static [f64]] = &[&[], &[0.5], &[0.0, 0.75]];
    const BU: &'static [f64] = &[2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0];
    const BE: &'static [f64] = &[];
}

/// Runge-Kutta 4th order integrator.
///
/// This integrator uses the classical 4th order Runge-Kutta method to numerically integrate the ODE system.
/// It calculates four intermediate values (k1, k2, k3, k4) to estimate the next step solution.
#[derive(Debug, Clone, Copy, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct RK4;

impl ButcherTableau for RK4 {
    const C: &'static [f64] = &[0.0, 0.5, 0.5, 1.0];
    const A: &'static [&'static [f64]] = &[&[], &[0.5], &[0.0, 0.5], &[0.0, 0.0, 1.0]];
    const BU: &'static [f64] = &[1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0];
    const BE: &'static [f64] = &[];
}

/// Ralston's 4th order integrator.
///
/// This fourth order method is known as minimum truncation error RK4.
#[derive(Debug, Clone, Copy, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct RALS4;

impl ButcherTableau for RALS4 {
    const C: &'static [f64] = &[0.0, 0.4, 0.45573725, 1.0];
    const A: &'static [&'static [f64]] = &[
        &[],
        &[0.4],
        &[0.29697761, 0.158575964],
        &[0.21810040, -3.050965616, 3.83286476],
    ];
    const BU: &'static [f64] = &[0.17476028, -0.55148066, 1.20553560, 0.17118478];
    const BE: &'static [f64] = &[];
}

/// Runge-Kutta 5th order integrator
///
/// This integrator uses the 5th order Runge-Kutta method to numerically integrate the ODE system.
#[derive(Debug, Clone, Copy, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct RK5;

impl ButcherTableau for RK5 {
    const C: &'static [f64] = &[0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
    const A: &'static [&'static [f64]] = &[
        &[],
        &[0.2],
        &[0.075, 0.225],
        &[44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0],
        &[
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
        ],
        &[
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
        ],
        &[
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
        ],
    ];
    const BU: &'static [f64] = &[
        5179.0 / 57600.0,
        0.0,
        7571.0 / 16695.0,
        393.0 / 640.0,
        -92097.0 / 339200.0,
        187.0 / 2100.0,
        1.0 / 40.0,
    ];
    const BE: &'static [f64] = &[];
}

// ┌─────────────────────────────────────────────────────────┐
//  Embedded Runge-Kutta
// └─────────────────────────────────────────────────────────┘
/// Bogacki-Shampine 3(2) method
///
/// This method is known as `ode23` in MATLAB.
///
/// # Member variables
///
/// - `tol`: The tolerance for the estimated error.
/// - `safety_factor`: The safety factor for the step size adjustment.
/// - `min_step_size`: The minimum step size.
/// - `max_step_size`: The maximum step size.
/// - `max_step_iter`: The maximum number of iterations per step.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct BS23 {
    pub tol: f64,
    pub safety_factor: f64,
    pub min_step_size: f64,
    pub max_step_size: f64,
    pub max_step_iter: usize,
}

impl Default for BS23 {
    fn default() -> Self {
        Self {
            tol: 1e-3,
            safety_factor: 0.9,
            min_step_size: 1e-6,
            max_step_size: 1e-1,
            max_step_iter: 100,
        }
    }
}

impl BS23 {
    pub fn new(
        tol: f64,
        safety_factor: f64,
        min_step_size: f64,
        max_step_size: f64,
        max_step_iter: usize,
    ) -> Self {
        Self {
            tol,
            safety_factor,
            min_step_size,
            max_step_size,
            max_step_iter,
        }
    }
}

impl ButcherTableau for BS23 {
    const C: &'static [f64] = &[0.0, 0.5, 0.75, 1.0];
    const A: &'static [&'static [f64]] = &[
        &[],
        &[0.5],
        &[0.0, 0.75],
        &[2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0],
    ];
    const BU: &'static [f64] = &[2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0];
    const BE: &'static [f64] = &[7.0 / 24.0, 0.25, 1.0 / 3.0, 0.125];

    fn tol(&self) -> f64 {
        self.tol
    }
    fn safety_factor(&self) -> f64 {
        self.safety_factor
    }
    fn min_step_size(&self) -> f64 {
        self.min_step_size
    }
    fn max_step_size(&self) -> f64 {
        self.max_step_size
    }
    fn max_step_iter(&self) -> usize {
        self.max_step_iter
    }
}

/// Runge-Kutta-Fehlberg 4/5th order integrator.
///
/// This integrator uses the Runge-Kutta-Fehlberg method, which is an adaptive step size integrator.
/// It calculates six intermediate values (k1, k2, k3, k4, k5, k6) to estimate the next step solution and the error.
/// The step size is automatically adjusted based on the estimated error to maintain the desired tolerance.
///
/// # Member variables
///
/// - `tol`: The tolerance for the estimated error.
/// - `safety_factor`: The safety factor for the step size adjustment.
/// - `min_step_size`: The minimum step size.
/// - `max_step_size`: The maximum step size.
/// - `max_step_iter`: The maximum number of iterations per step.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct RKF45 {
    pub tol: f64,
    pub safety_factor: f64,
    pub min_step_size: f64,
    pub max_step_size: f64,
    pub max_step_iter: usize,
}

impl Default for RKF45 {
    fn default() -> Self {
        Self {
            tol: 1e-6,
            safety_factor: 0.9,
            min_step_size: 1e-6,
            max_step_size: 1e-1,
            max_step_iter: 100,
        }
    }
}

impl RKF45 {
    pub fn new(
        tol: f64,
        safety_factor: f64,
        min_step_size: f64,
        max_step_size: f64,
        max_step_iter: usize,
    ) -> Self {
        Self {
            tol,
            safety_factor,
            min_step_size,
            max_step_size,
            max_step_iter,
        }
    }
}

impl ButcherTableau for RKF45 {
    const C: &'static [f64] = &[0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0];
    const A: &'static [&'static [f64]] = &[
        &[],
        &[0.25],
        &[3.0 / 32.0, 9.0 / 32.0],
        &[1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0],
        &[439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0],
        &[
            -8.0 / 27.0,
            2.0,
            -3544.0 / 2565.0,
            1859.0 / 4104.0,
            -11.0 / 40.0,
        ],
    ];
    const BU: &'static [f64] = &[
        16.0 / 135.0,
        0.0,
        6656.0 / 12825.0,
        28561.0 / 56430.0,
        -9.0 / 50.0,
        2.0 / 55.0,
    ];
    const BE: &'static [f64] = &[
        25.0 / 216.0,
        0.0,
        1408.0 / 2565.0,
        2197.0 / 4104.0,
        -1.0 / 5.0,
        0.0,
    ];

    fn tol(&self) -> f64 {
        self.tol
    }
    fn safety_factor(&self) -> f64 {
        self.safety_factor
    }
    fn min_step_size(&self) -> f64 {
        self.min_step_size
    }
    fn max_step_size(&self) -> f64 {
        self.max_step_size
    }
    fn max_step_iter(&self) -> usize {
        self.max_step_iter
    }
}

/// Dormand-Prince 5(4) method
///
/// This is an adaptive step size integrator based on a 5th order Runge-Kutta method with
/// 4th order embedded error estimation.
///
/// # Member variables
///
/// - `tol`: The tolerance for the estimated error.
/// - `safety_factor`: The safety factor for the step size adjustment.
/// - `min_step_size`: The minimum step size.
/// - `max_step_size`: The maximum step size.
/// - `max_step_iter`: The maximum number of iterations per step.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct DP45 {
    pub tol: f64,
    pub safety_factor: f64,
    pub min_step_size: f64,
    pub max_step_size: f64,
    pub max_step_iter: usize,
}

impl Default for DP45 {
    fn default() -> Self {
        Self {
            tol: 1e-6,
            safety_factor: 0.9,
            min_step_size: 1e-6,
            max_step_size: 1e-1,
            max_step_iter: 100,
        }
    }
}

impl DP45 {
    pub fn new(
        tol: f64,
        safety_factor: f64,
        min_step_size: f64,
        max_step_size: f64,
        max_step_iter: usize,
    ) -> Self {
        Self {
            tol,
            safety_factor,
            min_step_size,
            max_step_size,
            max_step_iter,
        }
    }
}

impl ButcherTableau for DP45 {
    const C: &'static [f64] = &[0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
    const A: &'static [&'static [f64]] = &[
        &[],
        &[0.2],
        &[0.075, 0.225],
        &[44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0],
        &[
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
        ],
        &[
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
        ],
        &[
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
        ],
    ];
    const BU: &'static [f64] = &[
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
        0.0,
    ];
    const BE: &'static [f64] = &[
        5179.0 / 57600.0,
        0.0,
        7571.0 / 16695.0,
        393.0 / 640.0,
        -92097.0 / 339200.0,
        187.0 / 2100.0,
        1.0 / 40.0,
    ];

    fn tol(&self) -> f64 {
        self.tol
    }
    fn safety_factor(&self) -> f64 {
        self.safety_factor
    }
    fn min_step_size(&self) -> f64 {
        self.min_step_size
    }
    fn max_step_size(&self) -> f64 {
        self.max_step_size
    }
    fn max_step_iter(&self) -> usize {
        self.max_step_iter
    }
}

/// Tsitouras 5(4) method
///
/// This is an adaptive step size integrator based on a 5th order Runge-Kutta method with
/// 4th order embedded error estimation, using the coefficients from Tsitouras (2011).
///
/// # Member variables
///
/// - `tol`: The tolerance for the estimated error.
/// - `safety_factor`: The safety factor for the step size adjustment.
/// - `min_step_size`: The minimum step size.
/// - `max_step_size`: The maximum step size.
/// - `max_step_iter`: The maximum number of iterations per step.
///
/// # References
///
/// - Ch. Tsitouras, Comput. Math. Appl. 62 (2011) 770-780.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct TSIT45 {
    pub tol: f64,
    pub safety_factor: f64,
    pub min_step_size: f64,
    pub max_step_size: f64,
    pub max_step_iter: usize,
}

impl Default for TSIT45 {
    fn default() -> Self {
        Self {
            tol: 1e-6,
            safety_factor: 0.9,
            min_step_size: 1e-6,
            max_step_size: 1e-1,
            max_step_iter: 100,
        }
    }
}

impl TSIT45 {
    pub fn new(
        tol: f64,
        safety_factor: f64,
        min_step_size: f64,
        max_step_size: f64,
        max_step_iter: usize,
    ) -> Self {
        Self {
            tol,
            safety_factor,
            min_step_size,
            max_step_size,
            max_step_iter,
        }
    }
}

impl ButcherTableau for TSIT45 {
    const C: &'static [f64] = &[0.0, 0.161, 0.327, 0.9, 0.9800255409045097, 1.0, 1.0];
    const A: &'static [&'static [f64]] = &[
        &[],
        &[Self::C[1]],
        &[Self::C[2] - 0.335480655492357, 0.335480655492357],
        &[
            Self::C[3] - (-6.359448489975075 + 4.362295432869581),
            -6.359448489975075,
            4.362295432869581,
        ],
        &[
            Self::C[4] - (-11.74888356406283 + 7.495539342889836 - 0.09249506636175525),
            -11.74888356406283,
            7.495539342889836,
            -0.09249506636175525,
        ],
        &[
            Self::C[5]
                - (-12.92096931784711 + 8.159367898576159
                    - 0.0715849732814010
                    - 0.02826905039406838),
            -12.92096931784711,
            8.159367898576159,
            -0.0715849732814010,
            -0.02826905039406838,
        ],
        &[
            Self::BU[0],
            Self::BU[1],
            Self::BU[2],
            Self::BU[3],
            Self::BU[4],
            Self::BU[5],
        ],
    ];
    const BU: &'static [f64] = &[
        0.09646076681806523,
        0.01,
        0.4798896504144996,
        1.379008574103742,
        -3.290069515436081,
        2.324710524099774,
        0.0,
    ];
    const BE: &'static [f64] = &[
        0.001780011052226,
        0.000816434459657,
        -0.007880878010262,
        0.144711007173263,
        -0.582357165452555,
        0.458082105929187,
        1.0 / 66.0,
    ];

    fn tol(&self) -> f64 {
        self.tol
    }
    fn safety_factor(&self) -> f64 {
        self.safety_factor
    }
    fn min_step_size(&self) -> f64 {
        self.min_step_size
    }
    fn max_step_size(&self) -> f64 {
        self.max_step_size
    }
    fn max_step_iter(&self) -> usize {
        self.max_step_iter
    }
}

/// Runge-Kutta-Fehlberg 7/8th order integrator.
///
/// This integrator uses the Runge-Kutta-Fehlberg 7(8) method, an adaptive step size integrator.
/// It evaluates f(x,y) thirteen times per step, using embedded 7th and 8th order
/// Runge-Kutta estimates to estimate the solution and the error.
/// The 7th order solution is propagated, and the difference between the 8th and 7th
/// order solutions is used for error estimation and step size control.
///
/// # Member variables
///
/// - `tol`: The tolerance for the estimated error.
/// - `safety_factor`: The safety factor for the step size adjustment.
/// - `min_step_size`: The minimum step size.
/// - `max_step_size`: The maximum step size.
/// - `max_step_iter`: The maximum number of iterations per step.
///
/// # References
/// - Meysam Mahooti (2025). [Runge-Kutta-Fehlberg (RKF78)](https://www.mathworks.com/matlabcentral/fileexchange/61130-runge-kutta-fehlberg-rkf78), MATLAB Central File Exchange.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct RKF78 {
    pub tol: f64,
    pub safety_factor: f64,
    pub min_step_size: f64,
    pub max_step_size: f64,
    pub max_step_iter: usize,
}

impl Default for RKF78 {
    fn default() -> Self {
        Self {
            tol: 1e-7, // Higher precision default for a higher-order method
            safety_factor: 0.9,
            min_step_size: 1e-10, // Smaller min step for higher order
            max_step_size: 1e-1,
            max_step_iter: 100,
        }
    }
}

impl RKF78 {
    pub fn new(
        tol: f64,
        safety_factor: f64,
        min_step_size: f64,
        max_step_size: f64,
        max_step_iter: usize,
    ) -> Self {
        Self {
            tol,
            safety_factor,
            min_step_size,
            max_step_size,
            max_step_iter,
        }
    }
}

impl ButcherTableau for RKF78 {
    const C: &'static [f64] = &[
        0.0,
        2.0 / 27.0,
        1.0 / 9.0,
        1.0 / 6.0,
        5.0 / 12.0,
        1.0 / 2.0,
        5.0 / 6.0,
        1.0 / 6.0,
        2.0 / 3.0,
        1.0 / 3.0,
        1.0,
        0.0, // k12 is evaluated at x[i]
        1.0, // k13 is evaluated at x[i]+h
    ];

    const A: &'static [&'static [f64]] = &[
        // k1
        &[],
        // k2
        &[2.0 / 27.0],
        // k3
        &[1.0 / 36.0, 3.0 / 36.0],
        // k4
        &[1.0 / 24.0, 0.0, 3.0 / 24.0],
        // k5
        &[20.0 / 48.0, 0.0, -75.0 / 48.0, 75.0 / 48.0],
        // k6
        &[1.0 / 20.0, 0.0, 0.0, 5.0 / 20.0, 4.0 / 20.0],
        // k7
        &[
            -25.0 / 108.0,
            0.0,
            0.0,
            125.0 / 108.0,
            -260.0 / 108.0,
            250.0 / 108.0,
        ],
        // k8
        &[
            31.0 / 300.0,
            0.0,
            0.0,
            0.0,
            61.0 / 225.0,
            -2.0 / 9.0,
            13.0 / 900.0,
        ],
        // k9
        &[
            2.0,
            0.0,
            0.0,
            -53.0 / 6.0,
            704.0 / 45.0,
            -107.0 / 9.0,
            67.0 / 90.0,
            3.0,
        ],
        // k10
        &[
            -91.0 / 108.0,
            0.0,
            0.0,
            23.0 / 108.0,
            -976.0 / 135.0,
            311.0 / 54.0,
            -19.0 / 60.0,
            17.0 / 6.0,
            -1.0 / 12.0,
        ],
        // k11
        &[
            2383.0 / 4100.0,
            0.0,
            0.0,
            -341.0 / 164.0,
            4496.0 / 1025.0,
            -301.0 / 82.0,
            2133.0 / 4100.0,
            45.0 / 82.0,
            45.0 / 164.0,
            18.0 / 41.0,
        ],
        // k12
        &[
            3.0 / 205.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -6.0 / 41.0,
            -3.0 / 205.0,
            -3.0 / 41.0,
            3.0 / 41.0,
            6.0 / 41.0,
            0.0,
        ],
        // k13
        &[
            -1777.0 / 4100.0,
            0.0,
            0.0,
            -341.0 / 164.0,
            4496.0 / 1025.0,
            -289.0 / 82.0,
            2193.0 / 4100.0,
            51.0 / 82.0,
            33.0 / 164.0,
            12.0 / 41.0,
            0.0,
            1.0,
        ],
    ];

    // Coefficients for the 7th order solution (propagated solution)
    // BU_i = BE_i (8th order) - ErrorCoeff_i
    // ErrorCoeff_i = [-41/840, 0, ..., 0, -41/840 (for k11), 41/840 (for k12), 41/840 (for k13)]
    const BU: &'static [f64] = &[
        41.0 / 420.0, // 41/840 - (-41/840)
        0.0,
        0.0,
        0.0,
        0.0,
        34.0 / 105.0,
        9.0 / 35.0,
        9.0 / 35.0,
        9.0 / 280.0,
        9.0 / 280.0,
        41.0 / 420.0,  // 41/840 - (-41/840)
        -41.0 / 840.0, // 0.0 - (41/840)
        -41.0 / 840.0, // 0.0 - (41/840)
    ];

    // Coefficients for the 8th order solution (used for error estimation)
    // These are from the y[i+1] formula in the MATLAB description
    const BE: &'static [f64] = &[
        41.0 / 840.0,
        0.0,
        0.0,
        0.0,
        0.0,
        34.0 / 105.0,
        9.0 / 35.0,
        9.0 / 35.0,
        9.0 / 280.0,
        9.0 / 280.0,
        41.0 / 840.0,
        0.0,
        0.0,
    ];

    fn tol(&self) -> f64 {
        self.tol
    }
    fn safety_factor(&self) -> f64 {
        self.safety_factor
    }
    fn min_step_size(&self) -> f64 {
        self.min_step_size
    }
    fn max_step_size(&self) -> f64 {
        self.max_step_size
    }
    fn max_step_iter(&self) -> usize {
        self.max_step_iter
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Gauss-Legendre 4th order
// └─────────────────────────────────────────────────────────┘

// Correct coefficients for 4th-order Gauss-Legendre method
const SQRT3: f64 = 1.7320508075688772;
const C1: f64 = 0.5 - SQRT3 / 6.0;
const C2: f64 = 0.5 + SQRT3 / 6.0;
const A11: f64 = 0.25;
const A12: f64 = 0.25 - SQRT3 / 6.0;
const A21: f64 = 0.25 + SQRT3 / 6.0;
const A22: f64 = 0.25;
const B1: f64 = 0.5;
const B2: f64 = 0.5;

/// Enum for implicit solvers.
///
/// This enum defines the available implicit solvers for the Gauss-Legendre 4th order integrator.
/// Currently, there are two options: fixed-point iteration and Broyden's method.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub enum ImplicitSolver {
    FixedPoint,
    Broyden,
    //TrustRegion(f64, f64),
}

/// Gauss-Legendre 4th order integrator.
///
/// This integrator uses the 4th order Gauss-Legendre Runge-Kutta method, which is an implicit integrator.
/// It requires solving a system of nonlinear equations at each step, which is done using the specified implicit solver (e.g., fixed-point iteration).
/// The Gauss-Legendre method has better stability properties compared to explicit methods, especially for stiff ODEs.
///
/// # Member variables
///
/// - `solver`: The implicit solver to use.
/// - `tol`: The tolerance for the implicit solver.
/// - `max_step_iter`: The maximum number of iterations for the implicit solver per step.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct GL4 {
    pub solver: ImplicitSolver,
    pub tol: f64,
    pub max_step_iter: usize,
}

impl Default for GL4 {
    fn default() -> Self {
        GL4 {
            solver: ImplicitSolver::FixedPoint,
            tol: 1e-8,
            max_step_iter: 100,
        }
    }
}

impl GL4 {
    pub fn new(solver: ImplicitSolver, tol: f64, max_step_iter: usize) -> Self {
        GL4 {
            solver,
            tol,
            max_step_iter,
        }
    }
}

impl ODEIntegrator for GL4 {
    #[allow(non_snake_case)]
    #[inline]
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64> {
        let n = y.len();
        //let sqrt3 = 3.0_f64.sqrt();
        //let c = 0.5 * (3.0 - sqrt3) / 6.0;
        //let d = 0.5 * (3.0 + sqrt3) / 6.0;
        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];

        // Initial guess for k1, k2.
        problem.rhs(t, y, &mut k1)?;
        k2.copy_from_slice(&k1);

        match self.solver {
            ImplicitSolver::FixedPoint => {
                // Fixed-point iteration
                let mut y1 = vec![0.0; n];
                let mut y2 = vec![0.0; n];

                for _ in 0..self.max_step_iter {
                    let k1_old = k1.clone();
                    let k2_old = k2.clone();

                    for i in 0..n {
                        y1[i] = y[i] + dt * (A11 * k1[i] + A12 * k2[i]);
                        y2[i] = y[i] + dt * (A11 * k1[i] + A12 * k2[i]);
                    }

                    // Compute new k1 and k2
                    problem.rhs(t + C1 * dt, &y1, &mut k1)?;
                    problem.rhs(t + C2 * dt, &y2, &mut k2)?;

                    // Check for convergence
                    let mut max_diff = 0f64;
                    for i in 0..n {
                        max_diff = max_diff.max((k1[i] - k1_old[i]).abs());
                        max_diff = max_diff.max((k2[i] - k2_old[i]).abs());
                    }

                    if max_diff < self.tol {
                        break;
                    }
                }
            }
            ImplicitSolver::Broyden => {
                let m = 2 * n;
                let mut U = vec![0.0; m];
                U[..n].copy_from_slice(&k1);
                U[n..].copy_from_slice(&k2);

                // F_vec = F(U)
                let mut F_vec = vec![0.0; m];
                compute_F(problem, t, y, dt, &U, &mut F_vec)?;

                // Initialize inverse Jacobian matrix
                let mut J_inv = eye(m);

                // Repeat Broyden's method
                for _ in 0..self.max_step_iter {
                    // delta = - J_inv * F_vec
                    let delta = (&J_inv * &F_vec).mul_scalar(-1.0);

                    // U <- U + delta
                    U.iter_mut()
                        .zip(delta.iter())
                        .for_each(|(u, d)| *u += *d);

                    let mut F_new = vec![0.0; m];
                    compute_F(problem, t, y, dt, &U, &mut F_new)?;

                    // If infinity norm of F_new is less than tol, break
                    if F_new.norm(Norm::LInf) < self.tol {
                        break;
                    }

                    // Residual: delta_F = F_new - F_vec
                    let delta_F = F_new.sub_vec(&F_vec);

                    // J_inv * delta_F
                    let J_inv_delta_F = &J_inv * &delta_F;

                    let denom = delta.dot(&J_inv_delta_F);
                    if denom.abs() < 1e-12 {
                        break;
                    }

                    // Broyden's "good" update for the inverse Jacobian
                    // J_inv <- J_inv + ((delta - J_inv * delta_F) * delta^T * J_inv) / denom
                    let delta_minus_J_inv_delta_F = delta.sub_vec(&J_inv_delta_F).to_col();
                    let delta_T_J_inv = &delta.to_row() * &J_inv;
                    let update = (delta_minus_J_inv_delta_F * delta_T_J_inv) / denom;
                    J_inv = J_inv + update;
                    F_vec = F_new;
                }

                k1.copy_from_slice(&U[..n]);
                k2.copy_from_slice(&U[n..]);
            }
        }

        for i in 0..n {
            y[i] += dt * (B1 * k1[i] + B2 * k2[i]);
        }

        Ok(dt)
    }
}

//// Helper function to compute the function F(U) for the implicit solver.
//// y1 = y + dt*(c*k1 + d*k2 - sqrt3/2*(k2-k1))
//// y2 = y + dt*(c*k1 + d*k2 + sqrt3/2*(k2-k1))
//#[allow(non_snake_case)]
//fn compute_F<P: ODEProblem>(
//    problem: &P,
//    t: f64,
//    y: &[f64],
//    dt: f64,
//    c: f64,
//    d: f64,
//    sqrt3: f64,
//    U: &[f64],
//    F: &mut [f64],
//) -> Result<()> {
//    let n = y.len();
//    let mut y1 = vec![0.0; n];
//    let mut y2 = vec![0.0; n];
//
//    for i in 0..n {
//        let k1 = U[i];
//        let k2 = U[n + i];
//        y1[i] = y[i] + dt * (c * k1 + d * k2 - sqrt3 * (k2 - k1) / 2.0);
//        y2[i] = y[i] + dt * (c * k1 + d * k2 + sqrt3 * (k2 - k1) / 2.0);
//    }
//
//    let mut f1 = vec![0.0; n];
//    let mut f2 = vec![0.0; n];
//    problem.rhs(t + c * dt, &y1, &mut f1)?;
//    problem.rhs(t + d * dt, &y2, &mut f2)?;
//
//    // F = [ k1 - f1, k2 - f2 ]
//    for i in 0..n {
//        F[i] = U[i] - f1[i];
//        F[n + i] = U[n + i] - f2[i];
//    }
//    Ok(())
//}
/// Helper function to compute the residual F(U) = U - f(y + dt*A*U)
#[allow(non_snake_case)]
fn compute_F<P: ODEProblem>(
    problem: &P,
    t: f64,
    y: &[f64],
    dt: f64,
    U: &[f64], // U is a concatenated vector [k1, k2]
    F: &mut [f64],
) -> Result<()> {
    let n = y.len();
    let (k1_slice, k2_slice) = U.split_at(n);

    let mut y1 = vec![0.0; n];
    let mut y2 = vec![0.0; n];

    for i in 0..n {
        y1[i] = y[i] + dt * (A11 * k1_slice[i] + A12 * k2_slice[i]);
        y2[i] = y[i] + dt * (A21 * k1_slice[i] + A22 * k2_slice[i]);
    }
    
    // F is an output parameter, its parts f1 and f2 are stored temporarily
    let (f1, f2) = F.split_at_mut(n);
    problem.rhs(t + C1 * dt, &y1, f1)?;
    problem.rhs(t + C2 * dt, &y2, f2)?;

    // Compute final residual F = [k1 - f1, k2 - f2]
    for i in 0..n {
        f1[i] = k1_slice[i] - f1[i];
        f2[i] = k2_slice[i] - f2[i];
    }
    Ok(())
}
