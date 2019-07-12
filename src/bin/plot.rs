extern crate peroxide;
use peroxide::*;

fn main() {
    let init_state = State::<f64>::new(0f64, c!(1), c!(0));

    let mut ode_solver = ExplicitODE::new(test_fn);

    ode_solver
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state)
        .set_step_size(0.01)
        .set_times(1000);

    let result = ode_solver.integrate();

    let x = result.col(0);
    let y = result.col(1);

    let mut plt = Plot2D::new();
    plt.set_domain(x)
        .insert_image(y)
        .set_title("Test Figure")
        .set_fig_size((10, 6))
        .set_dpi(300)
        .set_legends(vec!["RK4".to_owned()])
        .set_path("example_data/test_plot.png");

    plt.savefig();
}

fn test_fn(st: &mut State<f64>) {
    let x = st.param;
    let y = &st.value;
    let dy = &mut st.deriv;
    dy[0] = (5f64*x.powi(2) - y[0]) / (x + y[0]).exp();
}