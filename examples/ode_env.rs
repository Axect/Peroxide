use peroxide::fuga::*;

#[allow(unused_variables)]
fn main() -> Result<(), Box<dyn Error>> {
    let t = seq(0, 10, 1);
    let y = t
        .iter()
        .enumerate()
        .map(|(i, &t)| t.powi(5 - i as i32))
        .collect::<Vec<f64>>();

    let c = cubic_hermite_spline(&t, &y, Quadratic)?;

    let test_problem = Test { cs: c };
    let basic_ode_solver = BasicODESolver::new(RK4);
    let (t_vec, y_vec) = basic_ode_solver.solve(&test_problem, (0f64, 10f64), 0.01)?;
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
            .set_path("example_data/ode_env.png")
            .savefig()?;
    }

    Ok(())
}

struct Test {
    cs: CubicHermiteSpline,
}

impl ODEProblem for Test {
    fn initial_conditions(&self) -> Vec<f64> {
        vec![1f64]
    }

    fn rhs(&self, t: f64, _y: &[f64], dy: &mut [f64]) -> anyhow::Result<()> {
        Ok(dy[0] = self.cs.eval(t))
    }
}
