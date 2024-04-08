use thiserror::Error;
use ODEError::*;

use crate::{util::non_macro::eye, traits::math::LinearOp, traits::fp::FPVector};

pub trait ODEProblem {
    fn initial_conditions(&self) -> Vec<f64>;
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> Result<(), ODEError>;
}


pub trait ODEIntegrator {
    fn step<P: ODEProblem>(&self, problem: &P, t: f64, y: &mut [f64], dt: f64) -> Result<f64, ODEError>;
}

#[derive(Debug, Clone, Copy, Error)]
pub enum ODEError {
    #[error("constraint violation")]
    ConstraintViolation,
    #[error("reached maximum number of iterations per step")]
    ReachedMaxStepIter,
}

pub trait ODESolver {
    fn solve<P: ODEProblem>(&self, problem: &P, t_span: (f64, f64), dt: f64) -> Result<(Vec<f64>, Vec<Vec<f64>>), ODEError>;
}

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
        let mut y = problem.initial_conditions();
        let mut t_vec = vec![t];
        let mut y_vec = vec![y.clone()];

        while t < t_span.1 {
            let dt_step = self.integrator.step(problem, t, &mut y, dt)?;
            t += dt_step;
            t_vec.push(t);
            y_vec.push(y.clone());
        }

        Ok((t_vec, y_vec))
    }
}

// ┌─────────────────────────────────────────────────────────┐
//  Runge-Kutta
// └─────────────────────────────────────────────────────────┘
/// Runge-Kutta 4th order
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
/// Runge-Kutta-Fehlberg 4/5th order
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

// ┌─────────────────────────────────────────────────────────┐
//  Gauss-Legendre 4th order
// └─────────────────────────────────────────────────────────┘
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ImplicitSolver {
    FixedPoint,
    Broyden,
    TrustRegion(f64, f64),
}

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GL4 {
    solver: ImplicitSolver,
    max_iter: usize,
    tol: f64,
}

impl GL4 {
    pub fn new(solver: ImplicitSolver, max_iter: usize, tol: f64) -> Self {
        GL4 {
            solver,
            max_iter,
            tol,
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
                for _ in 0..self.max_iter {
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
            ImplicitSolver::Broyden => {
                // Broyden's method
                let mut b = eye(n);

                for _ in 0..self.max_iter {
                    problem.rhs(t + c * dt, &y1, &mut k1)?;
                    problem.rhs(t + d * dt, &y2, &mut k2)?;

                    let mut f1 = vec![0.0; n];
                    let mut f2 = vec![0.0; n];
                    for i in 0..n {
                        f1[i] = y1[i] - y[i] - dt * (c * k1[i] + d * k2[i] - sqrt3 * (k2[i] - k1[i]) / 2.0);
                        f2[i] = y2[i] - y[i] - dt * (c * k1[i] + d * k2[i] + sqrt3 * (k2[i] - k1[i]) / 2.0);
                    }

                    let neg_dy1 = b.apply(&f1);
                    let neg_dy2 = b.apply(&f2);

                    for i in 0..n {
                        y1[i] -= neg_dy1[i];
                        y2[i] -= neg_dy2[i];
                    }

                    let mut max_diff = 0f64;
                    for i in 0..n {
                        max_diff = max_diff.max(neg_dy1[i].abs()).max(neg_dy2[i].abs());
                    }

                    if max_diff < self.tol {
                        break;
                    }

                    let mut df1 = vec![0.0; n];
                    let mut df2 = vec![0.0; n];
                    for i in 0..n {
                        df1[i] = f1[i] + dt * (c * neg_dy1[i] + d * neg_dy2[i] - sqrt3 * (neg_dy2[i] - neg_dy1[i]) / 2.0);
                        df2[i] = f2[i] + dt * (c * neg_dy1[i] + d * neg_dy2[i] + sqrt3 * (neg_dy2[i] - neg_dy1[i]) / 2.0);
                    }

                    for i in 0..n {
                        for j in 0..n {
                            b[(i, j)] -= (df1[i] * neg_dy1[j] + df2[i] * neg_dy2[j]) / (neg_dy1[j].powi(2) + neg_dy2[j].powi(2));
                        }
                    }
                }
            }
            ImplicitSolver::TrustRegion(eta, rho) => {
                // Trust-region method
                let mut r = self.tol;
                let mut d1 = vec![0.0; n];
                let mut d2 = vec![0.0; n];
                let mut g1 = vec![0.0; n];
                let mut g2 = vec![0.0; n];
                let mut b1 = eye(n);
                let mut b2 = eye(n);

                for _ in 0..self.max_iter {
                    problem.rhs(t + c * dt, &y1, &mut k1)?;
                    problem.rhs(t + d * dt, &y2, &mut k2)?;

                    let mut f1 = vec![0.0; n];
                    let mut f2 = vec![0.0; n];
                    for i in 0..n {
                        f1[i] = y1[i] - y[i] - dt * (c * k1[i] + d * k2[i] - sqrt3 * (k2[i] - k1[i]) / 2.0);
                        f2[i] = y2[i] - y[i] - dt * (c * k1[i] + d * k2[i] + sqrt3 * (k2[i] - k1[i]) / 2.0);
                    }

                    let mut norm_f1 = 0.0;
                    let mut norm_f2 = 0.0;
                    for i in 0..n {
                        norm_f1 += f1[i] * f1[i];
                        norm_f2 += f2[i] * f2[i];
                        g1[i] = f1[i];
                        g2[i] = f2[i];
                    }
                    norm_f1 = norm_f1.sqrt();
                    norm_f2 = norm_f2.sqrt();

                    let mut inner_iter_count = 0;
                    while norm_f1 > self.tol || norm_f2 > self.tol {
                        let h1 = b1.apply(&g1).fmap(|x| -x);
                        let h2 = b2.apply(&g2).fmap(|x| -x);

                        let mut norm_h1 = 0.0;
                        let mut norm_h2 = 0.0;
                        for i in 0..n {
                            norm_h1 += h1[i] * h1[i];
                            norm_h2 += h2[i] * h2[i];
                        }
                        norm_h1 = norm_h1.sqrt();
                        norm_h2 = norm_h2.sqrt();

                        if norm_h1 <= r || norm_h2 <= r {
                            d1.copy_from_slice(&h1[..]);
                            d2.copy_from_slice(&h2[..]);
                        } else {
                            let sigma1 = (r / norm_h1).powf(2.0);
                            let sigma2 = (r / norm_h2).powf(2.0);
                            for i in 0..n {
                                d1[i] = sigma1 * h1[i];
                                d2[i] = sigma2 * h2[i];
                            }
                        }

                        let mut y1_new = vec![0.0; n];
                        let mut y2_new = vec![0.0; n];
                        for i in 0..n {
                            y1_new[i] = y1[i] + d1[i];
                            y2_new[i] = y2[i] + d2[i];
                        }

                        problem.rhs(t + c * dt, &y1_new, &mut k1)?;
                        problem.rhs(t + d * dt, &y2_new, &mut k2)?;

                        let mut f1_new = vec![0.0; n];
                        let mut f2_new = vec![0.0; n];
                        for i in 0..n {
                            f1_new[i] = y1_new[i] - y[i] - dt * (c * k1[i] + d * k2[i] - sqrt3 * (k2[i] - k1[i]) / 2.0);
                            f2_new[i] = y2_new[i] - y[i] - dt * (c * k1[i] + d * k2[i] + sqrt3 * (k2[i] - k1[i]) / 2.0);
                        }

                        let mut norm_f1_new = 0.0;
                        let mut norm_f2_new = 0.0;
                        for i in 0..n {
                            norm_f1_new += f1_new[i] * f1_new[i];
                            norm_f2_new += f2_new[i] * f2_new[i];
                        }
                        norm_f1_new = norm_f1_new.sqrt();
                        norm_f2_new = norm_f2_new.sqrt();

                        let rho1 = (norm_f1 * norm_f1 - norm_f1_new * norm_f1_new) / (norm_f1 * norm_f1 - (g1.iter().zip(&d1).map(|(g, d)| g * d).sum::<f64>()).abs());
                        let rho2 = (norm_f2 * norm_f2 - norm_f2_new * norm_f2_new) / (norm_f2 * norm_f2 - (g2.iter().zip(&d2).map(|(g, d)| g * d).sum::<f64>()).abs());

                        if rho1 > eta && rho2 > eta {
                            y1.copy_from_slice(&y1_new[..]);
                            y2.copy_from_slice(&y2_new[..]);
                            f1.copy_from_slice(&f1_new[..]);
                            f2.copy_from_slice(&f2_new[..]);
                            g1.copy_from_slice(&f1_new[..]);
                            g2.copy_from_slice(&f2_new[..]);
                            norm_f1 = norm_f1_new;
                            norm_f2 = norm_f2_new;
                            r = r.max(rho * r);
                        } else {
                            r *= rho;
                        }

                        if r < self.tol {
                            break;
                        }

                        for i in 0..n {
                            for j in 0..n {
                                b1[(i, j)] += (f1[i] - g1[i]) * d1[j] / (d1.iter().map(|d| d * d).sum::<f64>());
                                b2[(i, j)] += (f2[i] - g2[i]) * d2[j] / (d2.iter().map(|d| d * d).sum::<f64>());
                            }
                        }

                        inner_iter_count += 1;
                        if inner_iter_count >= self.max_iter {
                            break;
                        }
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
