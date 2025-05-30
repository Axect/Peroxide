use peroxide::fuga::*;

#[allow(unused_variables)]
fn main() -> Result<(), Box<dyn Error>> {
    let initial_conditions = vec![1f64];
    let rkf = RKF45::new(1e-4, 0.9, 1e-6, 1e-1, 100);
    let basic_ode_solver = BasicODESolver::new(rkf);
    let (t_vec, y_vec) = basic_ode_solver.solve(&Test, (0f64, 10f64), 0.01, &initial_conditions)?;
    let y_vec: Vec<f64> = y_vec.into_iter().flatten().collect();

    #[cfg(feature = "plot")]
    {
        let mut plt = Plot2D::new();
        plt.set_domain(t_vec)
            .insert_image(y_vec)
            .set_xlabel(r"$t$")
            .set_ylabel(r"$y$")
            .set_style(PlotStyle::Nature)
            .tight_layout()
            .set_dpi(600)
            .set_path("example_data/rkf45_test.png")
            .savefig()?;
    }

    Ok(())
}

struct Test;

impl ODEProblem for Test {
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) -> anyhow::Result<()> {
        Ok(dy[0] = (5f64 * t.powi(2) - y[0]) / (t + y[0]).exp())
    }
}
