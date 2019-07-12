use std::collections::HashMap;
use BoundaryCondition::Dirichlet;
use ExMethod::{Euler, RK4};
use ::{MutFP, Matrix, FPVector};
use ODEOptions::{BoundCond, InitCond, Method, StepSize, StopCond, Times};
use {Dual, Real};
use ::{zeros, cat};

/// Explicit ODE Methods
///
/// * Euler : Euler 1st Order
/// * RK4 : Runge-Kutta 4th Order
#[derive(Debug, Copy, Clone, Hash, PartialOrd, PartialEq, Eq)]
pub enum ExMethod {
    Euler,
    RK4,
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

    pub fn to_dual(&self) -> State<Dual> {
        State {
            param: self.param.to_dual(),
            value: self
                .value
                .clone()
                .into_iter()
                .map(|x| x.to_dual())
                .collect::<Vec<Dual>>(),
            deriv: self
                .deriv
                .clone()
                .into_iter()
                .map(|x| x.to_dual())
                .collect::<Vec<Dual>>(),
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

pub type ExUpdater = fn(&mut State<f64>);
pub type ImUpdater = fn(&mut State<Dual>);

/// ODE solver
///
/// * `Records` : Type of container to contain results
/// * `Param` : Type of parameter
/// * `ODEMethod` : Explicit or Implicit
pub trait ODE {
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
}

#[derive(Clone)]
pub struct ExplicitODE {
    state: State<f64>,
    func: fn(&mut State<f64>),
    step_size: f64,
    method: ExMethod,
    init_cond: State<f64>,
    bound_cond1: (State<f64>, BoundaryCondition),
    bound_cond2: (State<f64>, BoundaryCondition),
    stop_cond: fn(&Self) -> bool,
    times: usize,
    options: HashMap<ODEOptions, bool>,
}

impl ExplicitODE {
    pub fn new(f: ExUpdater) -> Self {
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
        }
    }

    pub fn get_state(&self) -> &State<f64> {
        &self.state
    }
}

impl ODE for ExplicitODE {
    type Records = Matrix;
    type Param = f64;
    type ODEMethod = ExMethod;

    fn mut_update(&mut self) {
        match self.method {
            Euler => {
                // Set Derivative from state
                (self.func)(&mut self.state);
                let dt = self.step_size;
                self.state
                    .value
                    .mut_zip_with(|x, y| x + y * dt, &self.state.deriv);
                self.state.param += dt;
            }
            RK4 => {
                let h = self.step_size;
                let h2 = h / 2f64;

                // Set Derivative from state
                let yn = self.state.value.clone();
                (self.func)(&mut self.state);

                let k1 = self.state.deriv.clone();
                let k1_add = k1.fmap(|x| x * h2);
                self.state.param += h2;
                self.state.value.mut_zip_with(|x, y| x + y, &k1_add);
                (self.func)(&mut self.state);

                let k2 = self.state.deriv.clone();
                let k2_add = k2.zip_with(|x, y| h2*x - y, &k1_add);
                self.state.value.mut_zip_with(|x, y| x + y, &k2_add);
                (self.func)(&mut self.state);

                let k3 = self.state.deriv.clone();
                let k3_add = k3.zip_with(|x, y| h*x - y, &k2_add);
                self.state.param += h2;
                self.state.value.mut_zip_with(|x, y| x + y, &k3_add);
                (self.func)(&mut self.state);

                let k4 = self.state.deriv.clone();

                for i in 0 .. k1.len() {
                    self.state.value[i] = yn[i] + (k1[i] + 2f64*k2[i] + 2f64*k3[i] + k4[i]) * h / 6f64;
                }
            }
        }
    }

    fn integrate(&mut self) -> Self::Records {
        assert!(self.check_enough(), "Not enough fields!");

        let mut result = zeros(self.times + 1, self.state.value.len() + 1);

        result.subs_row(0, cat(self.state.param, self.state.value.clone()));

        for i in 1 .. self.times+1 {
            self.mut_update();
            result.subs_row(i, cat(self.state.param,self.state.value.clone()));
        }

        result
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
            },
            None => { return false; }
        }

        // Step size
        match self.options.get(&StepSize) {
            Some(x) => {
                if !*x {
                    return false;
                }
            },
            None => { return false; }
        }

        // Initial or Boundary
        match self.options.get(&InitCond) {
            None => { return false; }
            Some(x) => {
                if !*x {
                    match self.options.get(&BoundCond) {
                        None => { return false; }
                        Some(_) => ()
                    }
                }
            }
        }

        // Set Time?
        match self.options.get(&Times) {
            None => { return false; }
            Some(x) => {
                if !*x {
                    return false;
                }
            }
        }
        true
    }
}