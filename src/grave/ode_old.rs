use structure::matrix::*;
use structure::vector::*;
use std::convert::Into;
pub use self::ODEMethod::*;
use numerical::runge_kutta::rk4;
use numerical::bdf::bdf1;
use structure::dual::*;
use numerical::gauss_legendre::gl4;
use util::non_macro::zeros;
use util::non_macro::cat;
use numerical::runge_kutta::one_step_rk;
use numerical::bdf::one_step_bdf1;
use numerical::gauss_legendre::k_newton;

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

pub fn solve_with_condition<F, T, G>(f: F, init_value: Vec<f64>, param_range: (T, T), step: f64, method: ODEMethod, cond: G) -> Matrix
where F: Fn(Dual, Vec<Dual>) -> Vec<Dual> + Copy,
      T: Into<f64> + Copy,
      G: Fn(Vec<f64>) -> bool + Copy,
{
    let t_start = param_range.0.into();
    let t_end = param_range.1.into();
    let num = ((t_end - t_start) / step).round() as usize;

    match method {
        RK4 => {
            let g = |t: f64, xs: Vec<f64>| f(dual(t, 0.), merge_dual(xs.clone(), vec![0f64;xs.len()])).values();
            let mut xs = init_value.clone();
            let mut t = t_start;
            let mut result = zeros(num+1, xs.len() + 1);
            result.subs_row(0, cat(t, xs.clone()));
            let mut key = 0usize;
            for i in 1 .. num + 1 {
                if !cond(xs.clone()) {
                    key = i - 1;
                    break;
                } else {
                    xs = one_step_rk(t, xs.clone(), g, step);
                    t += step;
                    result.subs_row(i, cat(t, xs.clone()));
                }
            }
            result.take(key, Row)
        },
        BDF1(rtol) => {
            let mut t = t_start;
            let mut xs = init_value.clone();
            let mut records = zeros(num + 1, xs.len() + 1);
            let mut key = 0usize;
            for i in 0 .. (num+1) {
                if !cond(xs.clone()) {
                    key =i;
                    break;
                }
                records.subs_row(i, cat(t, xs.clone()));
                t += step;
                xs = one_step_bdf1(t, xs.clone(), f, step, rtol);
            }
            records.take(key, Row)
        },
        GL4(rtol) => {
            let mut t = t_start;
            let mut y_curr = init_value.clone();
            let mut records = zeros(num + 1, y_curr.len() + 1);
            records.subs_row(0, cat(t, y_curr.clone()));
            let mut key = 0usize;

            for i in 0 .. num {
                if !cond(y_curr.clone()) {
                    key = i;
                    break;
                }
                let (k1, k2) = k_newton(f, t, y_curr.clone(), step, rtol);
                y_curr = y_curr.add(&k1.fmap(|x| 0.5 * x * step).add(&k2.fmap(|x| 0.5 * x * step)));
                t += step;
                records.subs_row(i+1, cat(t, y_curr.clone()))
            }

            records.take(key, Row)
        }
    }
}