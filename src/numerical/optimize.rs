//! To optimize parametric model (non-linear regression)
//!
//! ## `Optimizer` structure
//!
//! ### Declaration
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::{Number, OptMethod, OptOption};
//! use std::collections::HashMap;
//!
//! pub struct Optimizer {
//!     domain: Vec<f64>,
//!     observed: Vec<f64>,
//!     func: fn(&Vec<f64>, Vec<Number>) -> Vec<Number>,
//!     param: Vec<Number>,
//!     max_iter: usize,
//!     error: f64,
//!     method: OptMethod,
//!     option: HashMap<OptOption, bool>,
//! }
//! ```
//!
//! ### Method (Should fill)
//!
//! * `new` : Declare new Optimizer. **Should be mutable**
//! * `set_init_param` : Input initial parameter
//! * `set_max_iter` : Set maximum number of iterations
//! * `set_method` : Set method to optimization
//!
//! ### Method (Optional)
//!
//! * `get_domain` : Get domain
//! * `get_error` : Root mean square error
//!
//! ### Method (Generate result)
//!
//! * `optimize` : Optimize
//!
//! ## Example
//!
//! * Optimize $y = x^n$ model with $y = x^2$ with gaussian noise.
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::*;
//!
//! fn main() {
//!     // To prepare noise
//!     let normal = Normal(0f64, 0.1f64);
//!     let normal2 = Normal(0f64, 100f64);
//!
//!     // Noise to domain
//!     let mut x = seq(0., 99., 1f64);
//!     x = zip_with(|a, b| (a + b).abs(), &x, &normal.sample(x.len()));
//!
//!     // Noise to image
//!     let mut y = x.fmap(|t| t.powi(2));
//!     y = zip_with(|a, b| a + b, &y, &normal2.sample(y.len()));
//!
//!     // Initial parameter
//!     let n_init = vec![1f64];
//!     let data = hstack!(x.clone(), y.clone());
//!
//!     // Optimizer setting
//!     let mut opt = Optimizer::new(data, quad);
//!     let p = opt.set_init_param(n_init)
//!         .set_max_iter(50)
//!         .set_method(LevenbergMarquardt)
//!         .optimize();
//!     p.print();                  // Optimized parameter
//!     opt.get_error().print();    // Optimized RMSE
//!
//!     // To prepare plotting
//!     let z = quad(&x, NumberVector::from_f64_vec(p)).to_f64_vec();
//!
//!     // Plot
//!     //#[cfg(feature = "plot")]
//!     //{
//!     //    let mut plt = Plot2D::new();
//!     //    plt.set_domain(x)
//!     //        .insert_image(y)    // plot data
//!     //        .insert_image(z)    // plot fit
//!     //        .set_legend(vec!["Data", "Fit"])
//!     //        .set_title("Test ($y=x^2$ with noise)")
//!     //        .set_path("example_data/lm_test.png")
//!     //        .set_marker(vec![Point, Line])
//!     //        .savefig().expect("Can't draw a plot");
//!     //}
//! }
//!
//! fn quad(x: &Vec<f64>, n: Vec<Number>) -> Vec<Number> {
//!     x.clone().into_iter()
//!         .map(|t| Number::from_f64(t))
//!         .map(|t| t.powf(n[0]))
//!         .collect()
//! }
//! ```
//!
//! ![LM test](https://raw.githubusercontent.com/Axect/Peroxide/master/example_data/lm_test.png)

pub use self::OptMethod::{GaussNewton, GradientDescent, LevenbergMarquardt};
use self::OptOption::{InitParam, MaxIter};
use numerical::utils::jacobian;
use operation::number::{Number, NumberVector};
use std::collections::HashMap;
use structure::matrix::{LinearAlgebra, LinearOps, Matrix};
use util::useful::max;

#[derive(Debug, Clone, Copy)]
pub enum OptMethod {
    GradientDescent(f64),
    GaussNewton,
    LevenbergMarquardt,
}

#[derive(Debug, Clone, Copy, PartialOrd, PartialEq, Eq, Hash)]
pub enum OptOption {
    InitParam,
    MaxIter,
}

/// Optimizer for optimization (non-linear regression)
///
/// # Methods
/// * Gradient Descent
/// * Gauss Newton
/// * Levenberg Marquardt
pub struct Optimizer<F> where F: Fn(&Vec<f64>, Vec<Number>) -> Vec<Number> {
    domain: Vec<f64>,
    observed: Vec<f64>,
    func: Box<F>,
    param: Vec<Number>,
    max_iter: usize,
    error: f64,
    method: OptMethod,
    option: HashMap<OptOption, bool>,
}

impl<F> Optimizer<F> where F: Fn(&Vec<f64>, Vec<Number>) -> Vec<Number> {
    pub fn new(data: Matrix, func: F) -> Self {
        let mut default_option: HashMap<OptOption, bool> = HashMap::new();
        default_option.insert(InitParam, false);
        default_option.insert(MaxIter, false);

        Optimizer {
            domain: data.col(0),
            observed: data.col(1),
            func: Box::new(func),
            param: vec![],
            max_iter: 0,
            error: 0f64,
            method: LevenbergMarquardt,
            option: default_option,
        }
    }

    pub fn get_domain(&self) -> Vec<f64> {
        self.domain.clone()
    }

    pub fn get_error(&self) -> f64 {
        self.error
    }

    pub fn set_init_param(&mut self, p: Vec<f64>) -> &mut Self {
        if let Some(x) = self.option.get_mut(&InitParam) {
            *x = true
        }

        self.param = NumberVector::from_f64_vec(p);
        self
    }

    pub fn set_max_iter(&mut self, n: usize) -> &mut Self {
        if let Some(x) = self.option.get_mut(&MaxIter) {
            *x = true
        }

        self.max_iter = n;
        self
    }

    pub fn set_method(&mut self, method: OptMethod) -> &mut Self {
        self.method = method;
        self
    }

    pub fn optimize(&mut self) -> Vec<f64> {
        // Receive initial data
        let (x_vec, y_vec) = (self.domain.clone(), self.observed.clone());
        let (p_init, max_iter) = (self.param.clone(), self.max_iter);
        let f = |p: Vec<Number>| (self.func)(&x_vec, p);

        // Take various form of initial data
        let p_init_vec = p_init.to_f64_vec();
        let y = y_vec.to_matrix();

        // Declare mutable values
        let mut p = p_init_vec.to_matrix();
        let mut j = jacobian(f, &p_init_vec);
        let mut y_hat = f(p_init.clone()).to_f64_vec().to_matrix();
        let mut jtj = &j.t() * &j;

        match self.method {
            GradientDescent(alpha) => {
                for _i in 0..max_iter {
                    let h = alpha * j.t() * (&y - &y_hat);
                    p = &p + &h;
                    j = jacobian(f, &p.data);
                    y_hat = f(NumberVector::from_f64_vec(p.data.clone()))
                        .to_f64_vec()
                        .to_matrix();
                }
            }

            GaussNewton => unimplemented!(),

            LevenbergMarquardt => {
                let mut chi2 = ((&y - &y_hat).t() * (&y - &y_hat))[(0, 0)];
                let mut nu = 2f64;
                let lambda_0 = 1e-3;
                let mut lambda = lambda_0 * max(jtj.diag());

                for _i in 0..max_iter {
                    let h: Matrix;

                    match (jtj.clone() + lambda * jtj.to_diag()).inv() {
                        Some(b) => h = b * j.t() * (&y - &y_hat),
                        None => break,
                    }

                    let p_temp = &p + &h;
                    let j_temp = jacobian(f, &p.data);
                    let y_hat_temp = f(NumberVector::from_f64_vec(p_temp.data.clone()))
                        .to_f64_vec()
                        .to_matrix();
                    let chi2_temp = ((&y - &y_hat_temp).t() * (&y - &y_hat_temp))[(0, 0)];

                    let rho = (chi2 - chi2_temp)
                        / (h.t() * (lambda * jtj.to_diag() * h.clone() + j.t() * (&y - &y_hat)))
                            [(0, 0)];

                    if rho > 0f64 {
                        p = p_temp;
                        j = j_temp;
                        jtj = &j.t() * &j;
                        y_hat = y_hat_temp;
                        chi2 = chi2_temp;
                        lambda *= max(vec![1f64 / 3f64, 1f64 - (2f64 * rho - 1f64).powi(3)]);
                        nu = 2f64;
                    } else {
                        lambda *= nu;
                        nu *= 2f64;
                    }
                }
            }
        }
        let error_temp = &y - &y_hat;
        self.error = ((error_temp.t() * error_temp)[(0, 0)] / y.row as f64).sqrt();
        p.data
    }
}
