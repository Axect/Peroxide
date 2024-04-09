//! # Ordinary Differential Equation (ODE) Solvers
//!
//! This module provides traits and structs for solving ordinary differential equations (ODEs).
//!
//! ## Overview
//!
//! - `ODEProblem`: Trait for defining an ODE problem.
//! - `ODEIntegrator`: Trait for ODE integrators.
//! - `ODESolver`: Trait for ODE solvers.
//!
//! ## Available integrators
//!
//! - **Explicit**
//!   - Runge-Kutta 4th order (RK4)
//!   - Runge-Kutta-Fehlberg 4/5th order (RKF45)
//!   - Dormand-Prince 4/5th order (DP45)
//! - **Implicit**
//!   - Gauss-Legendre 4th order (GL4)
//!
//! ## Available solvers
//!
//! - `BasicODESolver`: A basic ODE solver using a specified integrator.
//!
//! ## Example
//!
//! ```rust
//! use peroxide::fuga::*;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     let rkf = RKF45::new(1e-4, 0.9, 1e-6, 1e-1, 100);
//!     let basic_ode_solver = BasicODESolver::new(rkf);
//!     let (t_vec, y_vec) = basic_ode_solver.solve(
//!         &Test,
//!         (0f64, 10f64),
//!         0.01,
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
//! struct Test;
//!
//! impl ODEProblem for Test {
//!     fn initial_conditions(&self) -> Vec<f64> {
//!         vec![1f64]
//!     }
//!
//!     fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> Result<(), ODEError> {
//!         Ok(dy[0] = (5f64 * t.powi(2) - y[0]) / (t + y[0]).exp())
//!     }
//! }
//! ```

use thiserror::Error;
use ODEError::*;

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
///     fn initial_conditions(&self) -> Vec<f64> {
///         vec![1.0, 2.0]
///     }
///
///     fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> Result<(), ODEError> {
///         dy[0] = -0.5 * y[0];
///         dy[1] = y[0] - y[1];
///         Ok(())
///     }
/// }
/// ```
pub trait ODEProblem {
    fn initial_conditions(&self) -> Vec<f64>;
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> Result<(), ODEError>;
}


/// Trait for ODE integrators.
///
/// Implement this trait to define your own ODE integrator.
pub trait ODEIntegrator {
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64, ODEError>;
}


/// Enum for ODE errors.


#[derive(Debug, Clone, Copy, Error)]
pub enum ODEError {
    #[error("constraint violation")]
    ConstraintViolation,
    #[error("reached maximum number of iterations per step")]
    ReachedMaxStepIter,
}

/// Trait for ODE solvers.
///
/// Implement this trait to define your own ODE solver.
pub trait ODESolver {
    fn solve<P: ODEProblem>(&self, problem: &P, t_span: (f64, f64), dt: f64) -> Result<(Vec<f64>, Vec<Vec<f64>>), ODEError>;
}

/// A basic ODE solver using a specified integrator.
///
/// # Example
///
/// ```
/// use peroxide::fuga::*;
///
/// fn main() -> Result<(), Box<dyn Error>> {
///     let rkf = RKF45::new(1e-4, 0.9, 1e-6, 1e-1, 100);
///     let basic_ode_solver = BasicODESolver::new(rkf);
///     let (t_vec, y_vec) = basic_ode_solver.solve(
///         &Test,
///         (0f64, 10f64),
///         0.01,
///     )?;
///     let y_vec: Vec<f64> = y_vec.into_iter().flatten().collect();
///
///     Ok(())
/// }
///
/// struct Test;
///
/// impl ODEProblem for Test {
///     fn initial_conditions(&self) -> Vec<f64> {
///         vec![1f64]
///     }
///
///     fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> Result<(), ODEError> {
///         Ok(dy[0] = (5f64 * t.powi(2) - y[0]) / (t + y[0]).exp())
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
    fn solve<P: ODEProblem>(&self, problem: &P, t_span: (f64, f64), dt: f64) -> Result<(Vec<f64>, Vec<Vec<f64>>), ODEError> {
        let mut t = t_span.0;
        let mut dt = dt;
        let mut y = problem.initial_conditions();
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
//  Runge-Kutta
// └─────────────────────────────────────────────────────────┘
/// Runge-Kutta 4th order integrator.
///
/// This integrator uses the classical 4th order Runge-Kutta method to numerically integrate the ODE system.
/// It calculates four intermediate values (k1, k2, k3, k4) to estimate the next step solution.
#[derive(Debug, Clone, Copy, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RK4;

impl ODEIntegrator for RK4 {
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64, ODEError> {
        let n = y.len();
        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];
        let mut k3 = vec![0.0; n];
        let mut k4 = vec![0.0; n];

        problem.rhs(t, y, &mut k1)?;
        k1.iter_mut().for_each(|k| *k *= dt);

        let mut y_temp = y.to_vec();
        y_temp.iter_mut().zip(&k1).for_each(|(y, k)| *y += 0.5 * k);
        problem.rhs(t + 0.5 * dt, &y_temp, &mut k2)?;
        k2.iter_mut().for_each(|k| *k *= dt);

        y_temp.iter_mut().zip(&k2).for_each(|(y, k)| *y = y.to_owned() + 0.5 * k);
        problem.rhs(t + 0.5 * dt, &y_temp, &mut k3)?;
        k3.iter_mut().for_each(|k| *k *= dt);

        y_temp.iter_mut().zip(&k3).for_each(|(y, k)| *y = y.to_owned() + k);
        problem.rhs(t + dt, &y_temp, &mut k4)?;
        k4.iter_mut().for_each(|k| *k *= dt);

        y.iter_mut().zip(&k1).zip(&k2).zip(&k3).zip(&k4)
            .for_each(|((((y, k1), k2), k3), k4)| {
                *y += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            });

        Ok(dt)
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Runge-Kutta-Fehlberg
// └─────────────────────────────────────────────────────────┘
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
pub struct RKF45 {
    tol: f64,
    safety_factor: f64,
    min_step_size: f64,
    max_step_size: f64,
    max_step_iter: usize,
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
    pub fn new(tol: f64, safety_factor: f64, min_step_size: f64, max_step_size: f64, max_step_iter: usize) -> Self {
        Self {
            tol,
            safety_factor,
            min_step_size,
            max_step_size,
            max_step_iter,
        }
    }

    pub fn set_tol(&mut self, tol: f64) {
        self.tol = tol;
    }

    pub fn set_safety_factor(&mut self, safety_factor: f64) {
        self.safety_factor = safety_factor;
    }

    pub fn set_min_step_size(&mut self, min_step_size: f64) {
        self.min_step_size = min_step_size;
    }

    pub fn set_max_step_size(&mut self, max_step_size: f64) {
        self.max_step_size = max_step_size;
    }

    pub fn set_max_step_iter(&mut self, max_step_iter: usize) {
        self.max_step_iter = max_step_iter;
    }

    pub fn get_tol(&self) -> f64 {
        self.tol
    }

    pub fn get_safety_factor(&self) -> f64 {
        self.safety_factor
    }

    pub fn get_min_step_size(&self) -> f64 {
        self.min_step_size
    }

    pub fn get_max_step_size(&self) -> f64 {
        self.max_step_size
    }

    pub fn get_max_step_iter(&self) -> usize {
        self.max_step_iter
    }
}

impl ODEIntegrator for RKF45 {
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64, ODEError> {
        let mut iter_count = 0usize;
        let mut dt = dt;
        let n = y.len();

        loop {
            let mut k1 = vec![0.0; n];
            let mut k2 = vec![0.0; n];
            let mut k3 = vec![0.0; n];
            let mut k4 = vec![0.0; n];
            let mut k5 = vec![0.0; n];
            let mut k6 = vec![0.0; n];

            problem.rhs(t, y, &mut k1)?;
            let mut y_temp = y.to_vec();
            for i in 0..n {
                y_temp[i] = y[i] + dt * (1.0 / 4.0) * k1[i];
            }

            problem.rhs(t + dt * (1.0 / 4.0), &y_temp, &mut k2)?;
            for i in 0..n {
                y_temp[i] = y[i] + dt * (3.0 / 32.0 * k1[i] + 9.0 / 32.0 * k2[i]);
            }

            problem.rhs(t + dt * (3.0 / 8.0), &y_temp, &mut k3)?;
            for i in 0..n {
                y_temp[i] = y[i] + dt * (1932.0 / 2197.0 * k1[i] - 7200.0 / 2197.0 * k2[i] + 7296.0 / 2197.0 * k3[i]);
            }

            problem.rhs(t + dt * (12.0 / 13.0), &y_temp, &mut k4)?;
            for i in 0..n {
                y_temp[i] = y[i] + dt * (439.0 / 216.0 * k1[i] - 8.0 * k2[i] + 3680.0 / 513.0 * k3[i] - 845.0 / 4104.0 * k4[i]);
            }

            problem.rhs(t + dt, &y_temp, &mut k5)?;
            for i in 0..n {
                y_temp[i] = y[i] + dt * (-8.0 / 27.0 * k1[i] + 2.0 * k2[i] - 3544.0 / 2565.0 * k3[i] + 1859.0 / 4104.0 * k4[i] - 11.0 / 40.0 * k5[i]);
            }

            problem.rhs(t + dt * 0.5, &y_temp, &mut k6)?;

            let mut error = 0.0;
            for i in 0..n {
                let y_err = dt * (1.0 / 360.0 * k1[i] - 128.0 / 4275.0 * k3[i] - 2197.0 / 75240.0 * k4[i] + 1.0 / 50.0 * k5[i] + 2.0 / 55.0 * k6[i]);
                error += y_err * y_err;
            }

            error = error.sqrt() / (n as f64).sqrt();

            let new_dt = self.safety_factor * dt * (self.tol * dt / error).powf(0.25);
            let new_dt = new_dt.max(self.min_step_size).min(self.max_step_size);

            if error <= self.tol {
                for i in 0..n {
                    y[i] += dt * (16.0 / 135.0 * k1[i] + 6656.0 / 12825.0 * k3[i] + 28561.0 / 56430.0 * k4[i] - 9.0 / 50.0 * k5[i] + 2.0 / 55.0 * k6[i]);
                }
                return Ok(new_dt);
            } else {
                iter_count += 1;
                if iter_count >= self.max_step_iter {
                    return Err(ReachedMaxStepIter);
                }
                dt = new_dt;
            }
        }
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
pub struct DP45 {
    tol: f64,
    safety_factor: f64,
    min_step_size: f64,
    max_step_size: f64,
    max_step_iter: usize,
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
    pub fn new(tol: f64, safety_factor: f64, min_step_size: f64, max_step_size: f64, max_step_iter: usize) -> Self {
        Self {
            tol,
            safety_factor,
            min_step_size,
            max_step_size,
            max_step_iter,
        }
    }
}

impl ODEIntegrator for DP45 {
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64, ODEError> {
        let mut iter_count = 0usize;
        let mut dt = dt;
        let n = y.len();

        loop {
            let mut k1 = vec![0.0; n];
            let mut k2 = vec![0.0; n];
            let mut k3 = vec![0.0; n];
            let mut k4 = vec![0.0; n];
            let mut k5 = vec![0.0; n];
            let mut k6 = vec![0.0; n];
            let mut k7 = vec![0.0; n];

            problem.rhs(t, y, &mut k1)?;

            let mut y_temp = y.to_vec();
            for i in 0..n {
                y_temp[i] = y[i] + dt * (1.0 / 5.0) * k1[i];
            }
            problem.rhs(t + dt * (1.0 / 5.0), &y_temp, &mut k2)?;

            for i in 0..n {
                y_temp[i] = y[i] + dt * (3.0 / 40.0 * k1[i] + 9.0 / 40.0 * k2[i]);
            }
            problem.rhs(t + dt * (3.0 / 10.0), &y_temp, &mut k3)?;

            for i in 0..n {
                y_temp[i] = y[i] + dt * (44.0 / 45.0 * k1[i] - 56.0 / 15.0 * k2[i] + 32.0 / 9.0 * k3[i]);
            }
            problem.rhs(t + dt * (4.0 / 5.0), &y_temp, &mut k4)?;

            for i in 0..n {
                y_temp[i] = y[i] + dt * (19372.0 / 6561.0 * k1[i] - 25360.0 / 2187.0 * k2[i] + 64448.0 / 6561.0 * k3[i] - 212.0 / 729.0 * k4[i]);
            }
            problem.rhs(t + dt * (8.0 / 9.0), &y_temp, &mut k5)?;

            for i in 0..n {
                y_temp[i] = y[i] + dt * (9017.0 / 3168.0 * k1[i] - 355.0 / 33.0 * k2[i] + 46732.0 / 5247.0 * k3[i] + 49.0 / 176.0 * k4[i] - 5103.0 / 18656.0 * k5[i]);
            }
            problem.rhs(t + dt, &y_temp, &mut k6)?;

            for i in 0..n {
                y_temp[i] = y[i] + dt * (35.0 / 384.0 * k1[i] + 500.0 / 1113.0 * k3[i] + 125.0 / 192.0 * k4[i] - 2187.0 / 6784.0 * k5[i] + 11.0 / 84.0 * k6[i]);
            }
            problem.rhs(t + dt, &y_temp, &mut k7)?;

            let mut error = 0.0;
            for i in 0..n {
                let y_err = dt * (71.0 / 57600.0 * k1[i] - 71.0 / 16695.0 * k3[i] + 71.0 / 1920.0 * k4[i] - 17253.0 / 339200.0 * k5[i] + 22.0 / 525.0 * k6[i] - 1.0 / 40.0 * k7[i]);
                error += y_err * y_err;
            }
            error = (error / (n as f64)).sqrt();

            let new_dt = self.safety_factor * dt * (self.tol / error).powf(0.2);
            let new_dt = new_dt.max(self.min_step_size).min(self.max_step_size);

            if error <= self.tol {
                y.copy_from_slice(&y_temp);
                return Ok(new_dt);
            } else {
                iter_count += 1;
                if iter_count >= self.max_step_iter {
                    return Err(ReachedMaxStepIter);
                }
                dt = new_dt;
            }
        }
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Gauss-Legendre 4th order
// └─────────────────────────────────────────────────────────┘
/// Enum for implicit solvers.
///
/// This enum defines the available implicit solvers for the Gauss-Legendre 4th order integrator.
/// Currently, only the fixed-point iteration method is implemented.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ImplicitSolver {
    FixedPoint,
    //Broyden,
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
pub struct GL4 {
    solver: ImplicitSolver,
    tol: f64,
    max_step_iter: usize,
}

impl Default for GL4 {
    fn default() -> Self {
        GL4 {
            solver: ImplicitSolver::FixedPoint,
            tol: 1e-6,
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
    #[inline]
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64, ODEError> {
        let n = y.len();
        let sqrt3 = 3.0_f64.sqrt();
        let c = 0.5 * (3.0 - sqrt3) / 6.0;
        let d = 0.5 * (3.0 + sqrt3) / 6.0;

        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];

        let mut y1 = vec![0.0; n];
        let mut y2 = vec![0.0; n];

        match self.solver {
            ImplicitSolver::FixedPoint => {
                // Fixed-point iteration
                for _ in 0..self.max_step_iter {
                    problem.rhs(t + c * dt, &y1, &mut k1)?;
                    problem.rhs(t + d * dt, &y2, &mut k2)?;

                    let mut max_diff = 0f64;
                    for i in 0..n {
                        let y1_new = y[i] + dt * (c * k1[i] + d * k2[i] - sqrt3 * (k2[i] - k1[i]) / 2.0);
                        let y2_new = y[i] + dt * (c * k1[i] + d * k2[i] + sqrt3 * (k2[i] - k1[i]) / 2.0);
                        max_diff = max_diff.max((y1_new - y1[i]).abs()).max((y2_new - y2[i]).abs());
                        y1[i] = y1_new;
                        y2[i] = y2_new;
                    }

                    if max_diff < self.tol {
                        break;
                    }
                }
            }
        }

        for i in 0..n {
            y[i] += dt * 0.5 * (k1[i] + k2[i]);
        }

        Ok(dt)
    }
}
