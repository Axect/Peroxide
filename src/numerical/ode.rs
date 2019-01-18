use structure::matrix::*;
use structure::vector::*;
use std::ops::{Add, Sub, Mul, Div};
use std::convert::Into;
pub use self::ODEMethod::*;
use numerical::runge_kutta::rk4;
use numerical::bdf::bdf1;
use structure::dual::*;
use numerical::gauss_legendre::gl4;

#[derive(Debug, Copy, Clone)]
pub enum ODEMethod {
    RK4,
    BDF1(f64),
    GL4(f64),
}

/// ODE Solver
///
/// # Parameters
/// * `f = f(t, y)`
/// * `init_value = y_start`
/// * `param_range = (t_start, t_end)`
///
/// # Type
/// `solve: (F, Vec<f64>, (T, T), f64, ODEMethod) -> Matrix where Fn(Dual, Vec<Dual>) -> Vec<Dual> + Copy, T: Into<f64> + Copy`
///
/// # Methods
/// * `RK4`: Explicit Runge-Kutta 4th order
/// * `BDF1`: Backward Differentiation Formula 1st order (Backward Euler)
/// * `GL4`: Gauss-Legendre 4th order
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let init_val = c!(2, 1);
/// let result = solve(lotka_volterra, init_val, (0, 10), 1e-3, RK4);
/// result.print();
///
/// fn lotka_volterra(t: Dual, xs: Vec<Dual>) -> Vec<Dual> {
///     let a = 4.;
///     let c = 1.;
///
///     let x = xs[0];
///     let y = xs[1];
///
///     vec![
///         a * (x - x*y),
///         -c * (y - x*y)
///     ]
/// }
/// ```
pub fn solve<F, T>(f: F, init_value: Vec<f64>, param_range: (T, T), step: f64, method: ODEMethod) -> Matrix
    where F: Fn(Dual, Vec<Dual>) -> Vec<Dual> + Copy,
          T: Into<f64> + Copy
{
    let t_start = param_range.0.into();
    let t_end = param_range.1.into();
    let num = ((t_end - t_start) / step).round() as usize;

    match method {
        RK4 => {
            let g = |t: f64, xs: Vec<f64>| f(dual(t, 0.), merge_dual(xs.clone(), vec![0f64;xs.len()])).values();
            let result = rk4(
                t_start,
                init_value.into_iter().map(|p| p.into()).collect::<Vec<f64>>(),
                g,
                step,
                num,
            );
            result
        },
        BDF1(rtol) => {
            let result = bdf1(t_start, init_value, f, step, rtol, num);
            result
        },
        GL4(rtol) => {
            let result = gl4(f, t_start, init_value, step, rtol, num);
            result
        }
    }
}
