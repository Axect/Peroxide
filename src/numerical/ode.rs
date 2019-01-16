use structure::matrix::*;
use structure::vector::*;
use std::ops::{Add, Sub, Mul, Div};
pub use self::ODEMethod::*;

#[derive(Debug, Copy, Clone)]
pub enum ODEMethod {
    RK4,
    BDF1,
}

/// ODE Solver
///
/// # Parameters
/// * `f = f(t, y)`
/// * `init_value = y_start`
/// * `param_range = (t_start, t_end)`
pub fn solve<F, G>(f: F, init_value: Vec<G>, param_range: (f64, f64), step: f64, method: ODEMethod) -> Matrix
    where F: Fn(Vec<G>) -> Vec<G> + Copy,
          G: Add<f64> + Sub<f64> + Mul<f64> + Div<f64> + Copy,
{
    unimplemented!()
}
