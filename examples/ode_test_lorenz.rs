extern crate peroxide;
use peroxide::*;

fn main() {
    let mut ex_test = ExplicitODE::new(f);

    let init_state: State<f64> = State::new(0.0, vec![10.0, 1.0, 1.0], vec![0.0, 0.0, 0.0]);
    ex_test
        .set_initial_condition(init_state)
        .set_method(ExMethod::Euler)
        .set_step_size(0.01f64)
        .set_times(10000);

    let mut ex_test2 = ex_test.clone();

    let _results = ex_test.integrate();

    ex_test2.set_method(ExMethod::RK4);

    let _results2 = ex_test2.integrate();

    // let mut wt = SimpleWriter::new();
    // wt.set_path("example_data/lorenz.pickle")
    //     .insert_matrix(results)
    //     .insert_matrix(results2)
    //     .write_pickle();
}

fn f(st: &mut State<f64>) {
    let x = &st.value;
    let dx = &mut st.deriv;
    dx[0] = 10f64 * (x[1] - x[0]);
    dx[1] = 28f64 * x[0] - x[1] - x[0] * x[2];
    dx[2] = -8f64 / 3f64 * x[2] + x[0] * x[1];
}
