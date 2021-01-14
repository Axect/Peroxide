#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let mut im_test = ImplicitODE::new(test_fn);

    let init_state = State::<f64>::new(0f64, c!(1), c!(0));
    im_test
        .set_initial_condition(init_state)
        .set_method(ImMethod::GL4)
        .set_step_size(0.01f64)
        .set_rtol(1e-6)
        .set_times(1000);

    #[cfg(feature = "plot")]
    {
        let result = im_test.integrate();
        let x = result.col(0);
        let y = result.col(1);
        let mut plt = Plot2D::new();
        plt.set_domain(x)
            .insert_image(y)
            .set_title("Test Figure")
            .set_fig_size((10, 6))
            .set_dpi(300)
            .set_legend(vec!["GL4"])
            .set_path("example_data/gl4_plot.png");
        plt.savefig();
    }
}

fn test_fn(st: &mut State<AD>, _: &NoEnv) {
    let x = st.param;
    let y = &st.value;
    let dy = &mut st.deriv;
    dy[0] = (5f64 * x.powi(2) - y[0]) / (x + y[0]).exp();
}
