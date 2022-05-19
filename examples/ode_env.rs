#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let x = seq(0, 10, 1);
    x.print();
    let y = x
        .iter()
        .enumerate()
        .map(|(i, &t)| t.powi(5 - i as i32))
        .collect::<Vec<f64>>();

    let c = CubicSpline::from_nodes(&x, &y);

    let init_state = State::<f64>::new(0f64, c!(1), c!(0));

    let mut ode_solver = ExplicitODE::new(test_fn);

    ode_solver
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state)
        .set_step_size(0.01)
        .set_times(1000)
        .set_env(c);

    let result = ode_solver.integrate();
    result.print();
}

fn test_fn(st: &mut State<f64>, env: &CubicSpline) {
    let x = st.param;
    let dy = &mut st.deriv;
    dy[0] = env.eval(x);
}
