#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

//fn main() {
//    let init_state = State::<f64>::new(0f64, c!(1), c!(0));
//
//    let mut ode_solver = ExplicitODE::new(test_fn);
//
//    ode_solver
//        .set_method(ExMethod::RK4)
//        .set_initial_condition(init_state)
//        .set_step_size(0.01)
//        .set_times(1000);
//
//    let _result = ode_solver.integrate();
//
//    // let mut st = SimpleWriter::new();
//    // st.set_path("example_data/rk4_test.pickle")
//    //     .insert_matrix(result)
//    //     .write_pickle();
//}
//
//fn test_fn(st: &mut State<f64>, _: &NoEnv) {
//    let x = st.param;
//    let y = &st.value;
//    let dy = &mut st.deriv;
//    dy[0] = (5f64 * x.powi(2) - y[0]) / (x + y[0]).exp();
//}

fn main() -> Result<(), Box<dyn Error>> {
    let basic_ode_solver = BasicODESolver::new(RK4);
    let (t_vec, y_vec) = basic_ode_solver.solve(
        &Test,
        (0f64, 10f64),
        1e-3,
    )?;
    let y_vec = y_vec.into_iter().flatten().collect();

    let mut plt = Plot2D::new();
    plt
        .set_domain(t_vec)
        .insert_image(y_vec)
        .set_xlabel(r"$t$")
        .set_ylabel(r"$y$")
        .set_style(PlotStyle::Nature)
        .tight_layout()
        .set_dpi(600)
        .set_path("example_data/rk4_test.png")
        .savefig()?;

    Ok(())
}

struct Test;

impl ODEProblem for Test {
    fn initial_conditions(&self) -> Vec<f64> {
        vec![1f64]
    }

    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> Result<(), ODEError> {
        Ok(dy[0] = (5f64 * t.powi(2) - y[0]) / (t + y[0]).exp())
    }
}
