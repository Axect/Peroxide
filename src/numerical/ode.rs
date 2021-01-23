//! Solver for ordinary differential equations
//!
//! ## Introduce `ODE` Trait & Structure
//!
//! ### `ODE` Trait
//!
//! * `ODE` structures are divided by two kinds
//!     * `ExplicitODE`
//!     * `ImplicitODE`
//! * `ODE` trait is given as
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::{Real, State, BoundaryCondition, Environment};
//!
//!     pub trait ODE<E: Environment> {
//!         type Records;
//!         type Vector;
//!         type Param;
//!         type ODEMethod;
//!
//!         fn mut_update(&mut self);
//!         fn integrate(&mut self) -> Self::Records;
//!         fn set_initial_condition<T: Real>(&mut self, init: State<T>) -> &mut Self;
//!         fn set_boundary_condition<T: Real>(
//!             &mut self,
//!             bound1: (State<T>, BoundaryCondition),
//!             bound2: (State<T>, BoundaryCondition),
//!         ) -> &mut Self;
//!         fn set_step_size(&mut self, dt: f64) -> &mut Self;
//!         fn set_method(&mut self, method: Self::ODEMethod) -> &mut Self;
//!         fn set_stop_condition(&mut self, f: fn(&Self) -> bool) -> &mut Self;
//!         fn set_times(&mut self, n: usize) -> &mut Self;
//!         fn check_enough(&self) -> bool;
//!         fn set_env(&mut self, env: E) -> &mut Self;
//!     }
//!     ```
//!
//!     * `Records` : The type to save results of ODE. Usually `Matrix` is used.
//!     * `Vector` : Vector can be below things.
//!         * `Vec<f64>` : Used for `ExplicitODE`
//!         * `Vec<AD>` : Used for `ImplicitODE`
//!     * `Param` : Also it can be `f64` or `AD`
//!     * `ODEMethod` : Method for solving ODE
//!         * `ExMethod` : Explicit method
//!             * `Euler` : Euler first order
//!             * `RK4` : Runge Kutta 4th order
//!         * `ImMethod` : Implicit method **(to be implemented)**
//!             * `BDF` : Backward Euler 1st order
//!             * `GL4` : Gauss Legendre 4th order
//!     * `Environment` : External environment (CubicSpline, Vec<f64>, Matrix or Another external table)
//!
//!
//! ### `State<T>` structure
//!
//! * To use `ODE` trait, you should understand `State<T>` first.
//!
//!     ```rust
//!     extern crate peroxide;
//!     use peroxide::fuga::Real;
//!
//!     #[derive(Debug, Clone, Default)]
//!     pub struct State<T: Real> {
//!         pub param: T,
//!         pub value: Vec<T>,
//!         pub deriv: Vec<T>,
//!     }
//!     ```
//!
//!     * `T` can be `f64` or `AD`
//!     * `param` is parameter for ODE. Usually it is represented by time.
//!     * `value` is value of each node.
//!     * `deriv` is value of derivative of each node.
//!
//! For example,
//!
//! $$ \frac{dy_n}{dt} = f(t, y_n) $$
//!
//! * $t$ is `param`
//! * $y_n$ is `value`
//! * $f(t,y_n)$ is `deriv`
//!
//! Methods for `State<T>` are as follows.
//!
//! * `to_f64(&self) -> State<f64>`
//! * `to_ad(&self) -> State<AD>`
//! * `new(T, Vec<T>, Vec<T>) -> Self`
//!
//! ### `Environment`
//!
//! * `Environment` needs just `Default`
//! * To use custom `Environment`, just type follows : `impl Environment for CustomType {}`
//! * If you don't want to use `Environment`, then use `NoEnv`
//! * Implemented Data Types
//!     * `Vec<f64>`
//!     * `Polynomial`
//!     * `Matrix`
//!     * `CubicSpline`
//!     * `NoEnv`
//!
//! ```
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let x = seq(0, 10, 1);
//!     x.print();
//!     let y = x.iter().enumerate().map(|(i, &t)| t.powi(5-i as i32)).collect::<Vec<f64>>();
//!
//!     let c = CubicSpline::from_nodes(x, y);
//!
//!     let init_state = State::<f64>::new(0f64, c!(1), c!(0));
//!
//!     let mut ode_solver = ExplicitODE::new(test_fn);
//!
//!     ode_solver
//!         .set_method(ExMethod::RK4)
//!         .set_initial_condition(init_state)
//!         .set_step_size(0.01)
//!         .set_times(1000)
//!         .set_env(c);
//!
//!     let result = ode_solver.integrate();
//!     result.print();
//! }
//!
//! fn test_fn(st: &mut State<f64>, env: &CubicSpline) {
//!     let x = st.param;
//!     let dy = &mut st.deriv;
//!     dy[0] = env.eval(x);
//! }
//! ```
//!
//! ### `ExplicitODE` struct
//!
//! `ExplicitODE` is given as follow :
//!
//! ```rust
//! extern crate peroxide;
//! use std::collections::HashMap;
//! use peroxide::fuga::{State, ExMethod, BoundaryCondition, ODEOptions, Environment};
//!
//! #[derive(Clone)]
//! pub struct ExplicitODE<E: Environment> {
//!     state: State<f64>,
//!     func: fn(&mut State<f64>, &E),
//!     step_size: f64,
//!     method: ExMethod,
//!     init_cond: State<f64>,
//!     bound_cond1: (State<f64>, BoundaryCondition),
//!     bound_cond2: (State<f64>, BoundaryCondition),
//!     stop_cond: fn(&Self) -> bool,
//!     times: usize,
//!     to_use: HashMap<ODEOptions, bool>,
//!     env: E,
//! }
//! ```
//!
//! * `state` : Current param, value, derivative
//! * `func` : Function to update `state`
//! * `init_cond` : Initial condition
//! * `bound_cond1` : If boundary problem, then first boundary condition
//! * `bound_cond2` : second boundary condition
//! * `stop_cond` : Stop condition (stop before `times`)
//! * `times` : How many times do you want to update?
//! * `to_use` : Just check whether information is enough
//! * `env` : Environment
//!
//! ## Example
//!
//! ### Lorenz Butterfly
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // =========================================
//!     //  Declare ODE (Euler)
//!     // =========================================
//!     let mut ex_test = ExplicitODE::new(f);
//!
//!     let init_state: State<f64> = State::new(
//!         0.0,
//!         vec![10.0, 1.0, 1.0],
//!         vec![0.0, 0.0, 0.0],
//!     );
//!
//!     ex_test
//!         .set_initial_condition(init_state)
//!         .set_method(ExMethod::Euler)
//!         .set_step_size(0.01f64)
//!         .set_times(10000);
//!
//!     // =========================================
//!     //  Declare ODE (RK4)
//!     // =========================================
//!     let mut ex_test2 = ex_test.clone();
//!     ex_test2.set_method(ExMethod::RK4);
//!
//!     // =========================================
//!     //  Declare ODE (GL4)
//!     // =========================================
//!     let mut im_test = ImplicitODE::new(g);
//!
//!     let init_state: State<AD> = State::new(
//!         AD0(0.0),
//!         vec![AD0(10f64), AD0(1f64), AD0(1f64)],
//!         vec![AD0(0.0), AD0(0.0), AD0(0.0)],
//!     );
//!
//!     im_test
//!         .set_initial_condition(init_state)
//!         .set_method(ImMethod::GL4)
//!         .set_step_size(0.01f64)
//!         .set_times(10000);
//!
//!     // =========================================
//!     //  Save results
//!     // =========================================
//!     let results = ex_test.integrate();
//!     let results2 = ex_test2.integrate();
//!     let results3 = im_test.integrate();
//!
//!     // Extract data
//!     let mut df = DataFrame::new(vec![]);
//!     df.push("x_euler", Series::new(results.col(1)));
//!     df.push("z_euler", Series::new(results.col(3)));
//!     df.push("x_rk4", Series::new(results2.col(1)));
//!     df.push("z_rk4", Series::new(results2.col(3)));
//!     df.push("x_gl4", Series::new(results3.col(1)));
//!     df.push("z_gl4", Series::new(results3.col(3)));
//!
//!     # #[cfg(feature="nc")]
//!     # {
//!     // Write netcdf file (`nc` feature required)
//!     df.write_nc("example_data/lorenz.nc")
//!         .expect("Can't write lorenz.nc");
//!     # }
//! }
//!
//! fn f(st: &mut State<f64>, _: &NoEnv) {
//!     let x = &st.value;
//!     let dx = &mut st.deriv;
//!     dx[0] = 10f64 * (x[1] - x[0]);
//!     dx[1] = 28f64 * x[0] - x[1] - x[0] * x[2];
//!     dx[2] = -8f64/3f64 * x[2] + x[0] * x[1];
//! }
//!
//! fn g(st: &mut State<AD>, _: &NoEnv) {
//!     let x = &st.value;
//!     let dx = &mut st.deriv;
//!     dx[0] = 10f64 * (x[1] - x[0]);
//!     dx[1] = 28f64 * x[0] - x[1] - x[0] * x[2];
//!     dx[2] = -8f64/3f64 * x[2] + x[0] * x[1];
//! }
//! ```
//!
//! If plotting pickle data with python, then
//!
//! ![Lorenz with Euler](https://raw.githubusercontent.com/Axect/Peroxide/master/example_data/lorenz_euler.png)
//!
//! ![Lorenz with RK4](https://raw.githubusercontent.com/Axect/Peroxide/master/example_data/lorenz_rk4.png)
//!
//! ![Lorenz with GL4](https://raw.githubusercontent.com/Axect/Peroxide/master/example_data/lorenz_gl4.png)
//!
//! ### Simple 1D Runge-Kutta
//!
//! $$\begin{gathered} \frac{dy}{dx} = \frac{5x^2 - y}{e^{x+y}} \\\ y(0) = 1 \end{gathered}$$
//!
//! ```rust
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let init_state = State::<f64>::new(0f64, c!(1), c!(0));
//!
//!     let mut ode_solver = ExplicitODE::new(test_fn);
//!
//!     ode_solver
//!         .set_method(ExMethod::RK4)
//!         .set_initial_condition(init_state)
//!         .set_step_size(0.01)
//!         .set_times(1000);
//!
//!     let result = ode_solver.integrate();
//!
//!     // Plot or Extract..
//! }
//!
//! fn test_fn(st: &mut State<f64>, _: &NoEnv) {
//!     let x = st.param;
//!     let y = &st.value;
//!     let dy = &mut st.deriv;
//!     dy[0] = (5f64*x.powi(2) - y[0]) / (x + y[0]).exp();
//! }
//! ```

use self::BoundaryCondition::Dirichlet;
use self::ExMethod::{Euler, RK4};
use self::ImMethod::{BDF1, GL4};
use self::ODEOptions::{BoundCond, InitCond, Method, StepSize, StopCond, Times};
use crate::numerical::{spline::CubicSpline, utils::jacobian};
use crate::structure::{
    ad::{AD, ADVec, AD::*},
    matrix::{LinearAlgebra, Matrix},
    polynomial::Polynomial,
};
use crate::traits::{
    fp::{FPMatrix, FPVector},
    math::{Norm, Normed, Vector},
    mutable::MutFP,
    num::Real,
};
use crate::util::{
    non_macro::{cat, concat, eye, zeros},
    print::Printable,
};
use std::collections::HashMap;
//#[cfg(feature = "O3")]
//use {blas_daxpy, blas_daxpy_return};

/// Explicit ODE Methods
///
/// * Euler : Euler 1st Order
/// * RK4 : Runge-Kutta 4th Order
#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum ExMethod {
    Euler,
    RK4,
}

#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum ImMethod {
    BDF1,
    GL4,
}

/// Kinds of Boundary Conditions
///
/// * Dirichlet
/// * Neumann
#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum BoundaryCondition {
    Dirichlet,
    Neumann,
}

/// Options for ODE
///
/// * `InitCond` : Initial condition
/// * `BoundCond` : Boundary condition
/// * `Method` : methods of `ExMethod` or `ImMethod`
/// * `StopCond` : Stop condition
/// * `StepSize` : Step size
/// * `Times` : A number of times to integrate with specific step size
#[derive(Debug, Clone, Copy, Hash, PartialOrd, PartialEq, Eq)]
pub enum ODEOptions {
    InitCond,
    BoundCond,
    Method,
    StopCond,
    StepSize,
    Times,
}

pub trait Environment: Default {}

/// State for ODE
///
/// * `param` : Parameter of ODE (ex) time)
/// * `value` : Current value of ODE
/// * `deriv` : Current differential of values
#[derive(Debug, Clone, Default)]
pub struct State<T: Real> {
    pub param: T,
    pub value: Vec<T>,
    pub deriv: Vec<T>,
}

impl<T: Real> State<T> {
    pub fn to_f64(&self) -> State<f64> {
        State {
            param: self.param.to_f64(),
            value: self
                .value
                .clone()
                .into_iter()
                .map(|x| x.to_f64())
                .collect::<Vec<f64>>(),
            deriv: self
                .deriv
                .clone()
                .into_iter()
                .map(|x| x.to_f64())
                .collect::<Vec<f64>>(),
        }
    }

    pub fn to_ad(&self) -> State<AD> {
        State {
            param: self.param.to_ad(),
            value: self
                .value
                .clone()
                .into_iter()
                .map(|x| x.to_ad())
                .collect::<Vec<AD>>(),
            deriv: self
                .deriv
                .clone()
                .into_iter()
                .map(|x| x.to_ad())
                .collect::<Vec<AD>>(),
        }
    }

    pub fn new(param: T, state: Vec<T>, deriv: Vec<T>) -> Self {
        State {
            param,
            value: state,
            deriv,
        }
    }
}

/// ODE solver
///
/// * `Records` : Type of container to contain results
/// * `Param` : Type of parameter
/// * `ODEMethod` : Explicit or Implicit
pub trait ODE<E: Environment> {
    type Records;
    type Param;
    type ODEMethod;

    fn mut_update(&mut self);
    //fn mut_integrate(&mut self, rec: &mut Self::Records);
    fn integrate(&mut self) -> Self::Records;
    fn set_initial_condition<T: Real>(&mut self, init: State<T>) -> &mut Self;
    fn set_boundary_condition<T: Real>(
        &mut self,
        bound1: (State<T>, BoundaryCondition),
        bound2: (State<T>, BoundaryCondition),
    ) -> &mut Self;
    fn set_step_size(&mut self, dt: f64) -> &mut Self;
    fn set_method(&mut self, method: Self::ODEMethod) -> &mut Self;
    fn set_stop_condition(&mut self, f: fn(&Self) -> bool) -> &mut Self;
    fn set_times(&mut self, n: usize) -> &mut Self;
    fn check_enough(&self) -> bool;
    fn set_env(&mut self, env: E) -> &mut Self;
}

#[derive(Clone)]
pub struct ExplicitODE<E: Environment> {
    state: State<f64>,
    func: fn(&mut State<f64>, &E),
    step_size: f64,
    method: ExMethod,
    init_cond: State<f64>,
    bound_cond1: (State<f64>, BoundaryCondition),
    bound_cond2: (State<f64>, BoundaryCondition),
    stop_cond: fn(&Self) -> bool,
    times: usize,
    options: HashMap<ODEOptions, bool>,
    env: E,
}

impl<E: Environment> ExplicitODE<E> {
    pub fn new(f: fn(&mut State<f64>, &E)) -> Self {
        let mut default_to_use: HashMap<ODEOptions, bool> = HashMap::new();
        default_to_use.insert(InitCond, false);
        default_to_use.insert(StepSize, false);
        default_to_use.insert(BoundCond, false);
        default_to_use.insert(Method, false);
        default_to_use.insert(StopCond, false);
        default_to_use.insert(Times, false);

        ExplicitODE {
            state: Default::default(),
            func: f,
            step_size: 0.0,
            method: Euler,
            init_cond: Default::default(),
            bound_cond1: (Default::default(), Dirichlet),
            bound_cond2: (Default::default(), Dirichlet),
            stop_cond: |_x| false,
            times: 0,
            options: default_to_use,
            env: E::default(),
        }
    }

    pub fn get_state(&self) -> &State<f64> {
        &self.state
    }

    pub fn get_env(&self) -> &E {
        &self.env
    }
}

impl<E: Environment> ODE<E> for ExplicitODE<E> {
    type Records = Matrix;
    type Param = f64;
    type ODEMethod = ExMethod;

    fn mut_update(&mut self) {
        match self.method {
            Euler => {
                // Set Derivative from state
                (self.func)(&mut self.state, &self.env);
                let dt = self.step_size;

                //match () {
                //    #[cfg(feature = "oxidize")]
                //    () => {
                //        blas_daxpy(dt, &self.state.deriv, &mut self.state.value);
                //    }
                //    _ => {
                        self.state
                            .value
                            .mut_zip_with(|x, y| x + y * dt, &self.state.deriv);
                //    }
                //}
                self.state.param += dt;
            }
            RK4 => {
                let h = self.step_size;
                let h2 = h / 2f64;

                // Set Derivative from state
                let yn = self.state.value.clone();
                (self.func)(&mut self.state, &self.env);

                let k1 = self.state.deriv.clone();
                let k1_add = k1.mul_scalar(h2);
                self.state.param += h2;
                self.state.value.mut_zip_with(|x, y| x + y, &k1_add);
                (self.func)(&mut self.state, &self.env);

                let k2 = self.state.deriv.clone();
                let k2_add = k2.zip_with(|x, y| h2 * x - y, &k1_add);
                self.state.value.mut_zip_with(|x, y| x + y, &k2_add);
                (self.func)(&mut self.state, &self.env);

                let k3 = self.state.deriv.clone();
                let k3_add = k3.zip_with(|x, y| h * x - y, &k2_add);
                self.state.param += h2;
                self.state.value.mut_zip_with(|x, y| x + y, &k3_add);
                (self.func)(&mut self.state, &self.env);

                let k4 = self.state.deriv.clone();

                for i in 0..k1.len() {
                    self.state.value[i] =
                        yn[i] + (k1[i] + 2f64 * k2[i] + 2f64 * k3[i] + k4[i]) * h / 6f64;
                }
            }
        }
    }

    fn integrate(&mut self) -> Self::Records {
        assert!(self.check_enough(), "Not enough fields!");

        let mut result = zeros(self.times + 1, self.state.value.len() + 1);

        result.subs_row(0, &cat(self.state.param, &self.state.value));

        match self.options.get(&StopCond) {
            Some(stop) if *stop => {
                let mut key = 1usize;
                for i in 1..self.times + 1 {
                    self.mut_update();
                    result.subs_row(i, &cat(self.state.param, &self.state.value));
                    key += 1;
                    if (self.stop_cond)(&self) {
                        println!("Reach the stop condition!");
                        print!("Current values are: ");
                        cat(self.state.param, &self.state.value).print();
                        break;
                    }
                }
                return result.take_row(key);
            }
            _ => {
                for i in 1..self.times + 1 {
                    self.mut_update();
                    result.subs_row(i, &cat(self.state.param, &self.state.value));
                }
                return result;
            }
        }
    }

    fn set_initial_condition<T: Real>(&mut self, init: State<T>) -> &mut Self {
        if let Some(x) = self.options.get_mut(&InitCond) {
            *x = true
        }
        self.init_cond = init.to_f64();
        self.state = init.to_f64();
        self
    }

    fn set_boundary_condition<T: Real>(
        &mut self,
        bound1: (State<T>, BoundaryCondition),
        bound2: (State<T>, BoundaryCondition),
    ) -> &mut Self {
        if let Some(x) = self.options.get_mut(&BoundCond) {
            *x = true
        }
        self.bound_cond1 = (bound1.0.to_f64(), bound1.1);
        self.bound_cond2 = (bound2.0.to_f64(), bound2.1);
        self
    }

    fn set_step_size(&mut self, dt: f64) -> &mut Self {
        if let Some(x) = self.options.get_mut(&StepSize) {
            *x = true
        }
        self.step_size = dt;
        self
    }

    fn set_method(&mut self, method: Self::ODEMethod) -> &mut Self {
        if let Some(x) = self.options.get_mut(&Method) {
            *x = true
        }
        self.method = method;
        self
    }

    fn set_stop_condition(&mut self, f: fn(&Self) -> bool) -> &mut Self {
        if let Some(x) = self.options.get_mut(&StopCond) {
            *x = true
        }
        self.stop_cond = f;
        self
    }

    fn set_times(&mut self, n: usize) -> &mut Self {
        if let Some(x) = self.options.get_mut(&Times) {
            *x = true
        }
        self.times = n;
        self
    }

    fn check_enough(&self) -> bool {
        // Method
        match self.options.get(&Method) {
            Some(x) => {
                if !*x {
                    return false;
                }
            }
            None => {
                return false;
            }
        }

        // Step size
        match self.options.get(&StepSize) {
            Some(x) => {
                if !*x {
                    return false;
                }
            }
            None => {
                return false;
            }
        }

        // Initial or Boundary
        match self.options.get(&InitCond) {
            None => {
                return false;
            }
            Some(x) => {
                if !*x {
                    match self.options.get(&BoundCond) {
                        None => {
                            return false;
                        }
                        Some(_) => (),
                    }
                }
            }
        }

        // Set Time?
        match self.options.get(&Times) {
            None => {
                return false;
            }
            Some(x) => {
                if !*x {
                    return false;
                }
            }
        }
        true
    }

    fn set_env(&mut self, env: E) -> &mut Self {
        self.env = env;
        self
    }
}

#[derive(Clone)]
pub struct ImplicitODE<E: Environment> {
    state: State<AD>,
    func: fn(&mut State<AD>, &E),
    step_size: f64,
    rtol: f64,
    method: ImMethod,
    init_cond: State<f64>,
    bound_cond1: (State<f64>, BoundaryCondition),
    bound_cond2: (State<f64>, BoundaryCondition),
    stop_cond: fn(&Self) -> bool,
    times: usize,
    options: HashMap<ODEOptions, bool>,
    env: E,
}

impl<E: Environment> ImplicitODE<E> {
    pub fn new(f: fn(&mut State<AD>, &E)) -> Self {
        let mut default_to_use: HashMap<ODEOptions, bool> = HashMap::new();
        default_to_use.insert(InitCond, false);
        default_to_use.insert(StepSize, false);
        default_to_use.insert(BoundCond, false);
        default_to_use.insert(Method, false);
        default_to_use.insert(StopCond, false);
        default_to_use.insert(Times, false);

        ImplicitODE {
            state: State::new(AD0(0f64), vec![], vec![]),
            func: f,
            step_size: 0.0,
            rtol: 1e-6,
            method: GL4,
            init_cond: Default::default(),
            bound_cond1: (Default::default(), Dirichlet),
            bound_cond2: (Default::default(), Dirichlet),
            stop_cond: |_x| false,
            times: 0,
            options: default_to_use,
            env: E::default(),
        }
    }

    pub fn get_state(&self) -> &State<AD> {
        &self.state
    }

    pub fn set_rtol(&mut self, rtol: f64) -> &mut Self {
        self.rtol = rtol;
        self
    }

    pub fn get_env(&self) -> &E {
        &self.env
    }
}

/// Value of 3f64.sqrt()
const SQRT3: f64 = 1.7320508075688772;

/// Butcher tableau for Gauss_legendre 4th order
const GL4_TAB: [[f64; 3]; 2] = [
    [0.5 - SQRT3 / 6f64, 0.25, 0.25 - SQRT3 / 6f64],
    [0.5 + SQRT3 / 6f64, 0.25 + SQRT3 / 6f64, 0.25],
];

#[allow(non_snake_case)]
impl<E: Environment> ODE<E> for ImplicitODE<E> {
    type Records = Matrix;
    type Param = AD;
    type ODEMethod = ImMethod;
    fn mut_update(&mut self) {
        match self.method {
            BDF1 => unimplemented!(),
            GL4 => {
                let f = |t: AD, y: Vec<AD>| {
                    let mut st = State::new(t, y.clone(), y);
                    (self.func)(&mut st, &self.env);
                    st.deriv
                };

                let h = self.step_size;
                let t = self.state.param;
                let t1: AD = t + GL4_TAB[0][0] * h;
                let t2: AD = t + GL4_TAB[1][0] * h;
                let yn = &self.state.value;
                let n = yn.len();

                // 0. Initial Guess
                let k1_init: Vec<f64> = f(t, yn.clone()).to_f64_vec();
                let k2_init: Vec<f64> = f(t, yn.clone()).to_f64_vec();
                let mut k_curr: Vec<f64> = concat(&k1_init, &k2_init);
                let mut err = 1f64;

                // 1. Combine two functions to one function
                let g = |k: &Vec<AD>| -> Vec<AD> {
                    let k1 = k.take(n);
                    let k2 = k.skip(n);
                    concat(
                        &f(
                            t1,
                            yn.add_vec(
                                &k1.mul_scalar(AD::from(GL4_TAB[0][1] * h))
                                    .add_vec(&k2.mul_scalar(AD::from(GL4_TAB[0][2] * h))),
                            ),
                        ),
                        &f(
                            t2,
                            yn.add_vec(
                                &k1.mul_scalar(AD::from(GL4_TAB[1][1] * h))
                                    .add_vec(&k2.mul_scalar(AD::from(GL4_TAB[1][2] * h))),
                            ),
                        ),
                    )
                };

                // 2. Obtain Jacobian
                let I = eye(2 * n);

                let mut Dg = jacobian(Box::new(g), &k_curr);
                let mut DG = &I - &Dg;
                let mut DG_inv = DG.inv();
                let mut G = k_curr.sub_vec(&g(&k_curr.to_ad_vec()).to_f64_vec());
                let mut num_iter: usize = 0;

                // 3. Iteration
                while err >= self.rtol && num_iter <= 10 {
                    let DGG = &DG_inv * &G;
                    let k_prev = k_curr.clone();
                    k_curr.mut_zip_with(|x, y| x - y, &DGG);
                    Dg = jacobian(Box::new(g), &k_curr);
                    DG = &I - &Dg;
                    DG_inv = DG.inv();
                    G = k_curr.sub_vec(&g(&k_curr.to_ad_vec()).to_f64_vec());
                    err = k_curr.sub_vec(&k_prev).norm(Norm::L2);
                    num_iter += 1;
                }

                // 4. Find k1, k2
                let (k1, k2) = (k_curr.take(n), k_curr.skip(n));

                // Set Derivative from state
                (self.func)(&mut self.state, &self.env);

                // 5. Merge k1, k2
                let mut y_curr = self.state.value.to_f64_vec();
                y_curr = y_curr.add_vec(&k1.mul_scalar(0.5 * h).add_vec(&k2.mul_scalar(0.5 * h)));
                self.state.value = y_curr.to_ad_vec();
                self.state.param = self.state.param + h;
            }
        }
    }

    fn integrate(&mut self) -> Self::Records {
        assert!(self.check_enough(), "Not enough fields!");

        let mut result = zeros(self.times + 1, self.state.value.len() + 1);

        result.subs_row(
            0,
            &cat(self.state.param.to_f64(), &self.state.value.to_f64_vec()),
        );

        match self.options.get(&StopCond) {
            Some(stop) if *stop => {
                let mut key = 1usize;
                for i in 1..self.times + 1 {
                    self.mut_update();
                    result.subs_row(
                        i,
                        &cat(self.state.param.to_f64(), &self.state.value.to_f64_vec()),
                    );
                    key += 1;
                    if (self.stop_cond)(&self) {
                        println!("Reach the stop condition!");
                        print!("Current values are: ");
                        cat(self.state.param.to_f64(), &self.state.value.to_f64_vec()).print();
                        break;
                    }
                }
                return result.take_row(key);
            }
            _ => {
                for i in 1..self.times + 1 {
                    self.mut_update();
                    result.subs_row(
                        i,
                        &cat(self.state.param.to_f64(), &self.state.value.to_f64_vec()),
                    );
                }
                return result;
            }
        }
    }

    fn set_initial_condition<T: Real>(&mut self, init: State<T>) -> &mut Self {
        if let Some(x) = self.options.get_mut(&InitCond) {
            *x = true
        }
        self.init_cond = init.to_f64();
        self.state = init.to_ad();
        self
    }

    fn set_boundary_condition<T: Real>(
        &mut self,
        bound1: (State<T>, BoundaryCondition),
        bound2: (State<T>, BoundaryCondition),
    ) -> &mut Self {
        if let Some(x) = self.options.get_mut(&BoundCond) {
            *x = true
        }
        self.bound_cond1 = (bound1.0.to_f64(), bound1.1);
        self.bound_cond2 = (bound2.0.to_f64(), bound2.1);
        self
    }

    fn set_step_size(&mut self, dt: f64) -> &mut Self {
        if let Some(x) = self.options.get_mut(&StepSize) {
            *x = true
        }
        self.step_size = dt;
        self
    }

    fn set_method(&mut self, method: Self::ODEMethod) -> &mut Self {
        if let Some(x) = self.options.get_mut(&Method) {
            *x = true
        }
        self.method = method;
        self
    }

    fn set_stop_condition(&mut self, f: fn(&Self) -> bool) -> &mut Self {
        if let Some(x) = self.options.get_mut(&StopCond) {
            *x = true
        }
        self.stop_cond = f;
        self
    }

    fn set_times(&mut self, n: usize) -> &mut Self {
        if let Some(x) = self.options.get_mut(&Times) {
            *x = true
        }
        self.times = n;
        self
    }

    fn check_enough(&self) -> bool {
        // Method
        match self.options.get(&Method) {
            Some(x) => {
                if !*x {
                    return false;
                }
            }
            None => {
                return false;
            }
        }

        // Step size
        match self.options.get(&StepSize) {
            Some(x) => {
                if !*x {
                    return false;
                }
            }
            None => {
                return false;
            }
        }

        // Initial or Boundary
        match self.options.get(&InitCond) {
            None => {
                return false;
            }
            Some(x) => {
                if !*x {
                    match self.options.get(&BoundCond) {
                        None => {
                            return false;
                        }
                        Some(_) => (),
                    }
                }
            }
        }

        // Set Time?
        match self.options.get(&Times) {
            None => {
                return false;
            }
            Some(x) => {
                if !*x {
                    return false;
                }
            }
        }
        true
    }

    fn set_env(&mut self, env: E) -> &mut Self {
        self.env = env;
        self
    }
}

// =============================================================================
// Some Environments
// =============================================================================
#[derive(Debug, Copy, Clone, Default)]
pub struct NoEnv {}

impl Environment for NoEnv {}

impl Environment for CubicSpline {}

impl Environment for Matrix {}

impl Environment for Vec<f64> {}

impl Environment for Polynomial {}
