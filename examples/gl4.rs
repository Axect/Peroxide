extern crate peroxide;
use peroxide::fuga::*;

pub fn main() {
    let t = AD0(0f64);
    let y = AD1(5f64, 1f64);
    let init_state = State::<AD>::new(t, vec![y], vec![AD0(0f64)]);

    let mut ode_solver = ImplicitODE::new(f);

    ode_solver
        .set_method(ImMethod::GL4)
        .set_initial_condition(init_state)
        .set_step_size(1e-3)
        .set_times(10_000);

    let result = ode_solver.integrate();

    let mut df = DataFrame::new(vec![]);
    df.push("t", Series::new(result.col(0)));
    df.push("y", Series::new(result.col(1)));
    df.print();

}

fn f(st: &mut State<AD>, _: &NoEnv) {
    let y = st.value[0];
    let dy = &mut st.deriv;
    dy[0] = -0.3 * y;
}
