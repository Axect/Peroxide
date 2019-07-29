use std::collections::HashMap;
use ::{Matrix, Number, NumberVector, LinearAlgebra, LinearOps};
use ::{max, jacobian};
use ::OptOption::{InitParam, MaxIter};
pub use ::OptMethod::{GradientDescent, GaussNewton, LevenbergMarquardt};


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

pub struct Optimizer {
    domain: Vec<f64>,
    observed: Vec<f64>,
    func: fn(&Vec<f64>, Vec<Number>) -> Vec<Number>,
    param: Vec<Number>,
    max_iter: usize,
    method: OptMethod,
    option: HashMap<OptOption, bool>,
}

impl Optimizer {
    pub fn new(data: Matrix, func: fn(&Vec<f64>, Vec<Number>) -> Vec<Number>) -> Self {
        let mut default_option: HashMap<OptOption, bool> = HashMap::new();
        default_option.insert(InitParam, false);
        default_option.insert(MaxIter, false);

        Optimizer {
            domain: data.col(0),
            observed: data.col(1),
            func,
            param: vec![],
            max_iter: 0,
            method: LevenbergMarquardt,
            option: default_option,
        }
    }

    pub fn get_domain(&self) -> Vec<f64> {
        self.domain.clone()
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

    pub fn optimize(&self) -> Matrix {
        // Receive initial data
        let (x_vec, y_vec) = (self.domain.clone(), self.observed.clone());
        let (p_init, max_iter) = (self.param.clone(), self.max_iter);
        let f = |p: Vec<Number>| (self.func)(&x_vec, p);

        // Take various form of initial data
        let p_init_vec = p_init.to_f64_vec();
        let y = y_vec.to_matrix();

        // Declare mutable values
        let mut p = p_init_vec.to_matrix();
        let mut j = jacobian(f, p.data.clone());
        let mut y_hat = f(p_init.clone()).to_f64_vec().to_matrix();
        let mut jtj = &j.t() * &j;

        match self.method {
            GradientDescent(alpha) => {
                for _i in 0 .. max_iter {
                    let h = alpha * j.t() * (&y - &y_hat);
                    p = &p + &h;
                    j = jacobian(f, p.data.clone());
                    y_hat = f(NumberVector::from_f64_vec(p.data.clone())).to_f64_vec().to_matrix();
                }
            }

            GaussNewton => {
                unimplemented!()
            }

            LevenbergMarquardt => {
                let mut chi2 = ((&y - &y_hat).t() * (&y - &y_hat))[(0,0)];
                let mut nu = 2f64;
                let lambda_0 = 1e-2;
                let mut lambda = lambda_0 * max(jtj.diag());

                for _i in 0 .. max_iter {
                    let h: Matrix;

                    match (jtj.clone() + lambda * jtj.to_diag()).inv() {
                        Some(b) => h = b * j.t() * (&y - &y_hat),
                        None => break,
                    }

                    let p_temp = &p + &h;
                    let j_temp = jacobian(f, p.data.clone());
                    let y_hat_temp = f(NumberVector::from_f64_vec(p_temp.data.clone())).to_f64_vec().to_matrix();
                    let chi2_temp = ((&y - &y_hat_temp).t() * (&y - &y_hat_temp))[(0,0)];

                    let rho = (chi2 - chi2_temp) / (h.t() * (lambda * jtj.to_diag() * h.clone() + j.t() * (&y - &y_hat)))[(0,0)];

                    if rho > 0f64 {
                        p = p_temp;
                        j = j_temp;
                        jtj = &j.t() * &j;
                        y_hat = y_hat_temp;
                        chi2 = chi2_temp;
                        lambda *= max(vec![1f64/3f64, 1f64 - (2f64*rho - 1f64).powi(3)]);
                        nu = 2f64;
                    } else {
                        lambda *= nu;
                        nu *= 2f64;
                    }
                }
            }
        }
        p
    }
}

