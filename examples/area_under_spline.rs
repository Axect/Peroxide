use anyhow::Result;
use peroxide::fuga::*;

fn main() -> Result<()> {
    // 1. Define the B-spline
    let degree = 3;
    let knots = vec![0f64, 1f64, 2f64, 3f64, 4f64];
    let control_points = vec![
        vec![0f64, 0f64],
        vec![1f64, 2f64],
        vec![2f64, -1f64],
        vec![3f64, 1f64],
        vec![4f64, -1f64],
        vec![5f64, 2f64],
        vec![6f64, 0f64],
    ];

    let spline = BSpline::clamped(degree, knots, control_points.clone())?;
    let deriv_spline = spline.derivative()?;

    // 2. Define the integrand for ∫y dx = ∫ y(t) * (dx/dt) dt
    let integrand = |t: f64, s: &BSpline, ds: &BSpline| -> f64 {
        let p = s.eval(t);
        let dp = ds.eval(t);
        let y_t = p.1;
        let dx_dt = dp.0;
        y_t * dx_dt
    };

    // 3. Perform numerical integration using self-contained Gauss-Legendre quadrature
    let t_start = 0f64;
    let t_end = 4f64;

    //let area = gauss_legendre_integrate(
    //    |t| integrand(t, &spline, &deriv_spline),
    //    t_start,
    //    t_end,
    //    64,
    //);
    let area = integrate(
        |t| integrand(t, &spline, &deriv_spline),
        (t_start, t_end),
        G7K15R(1e-4, 20),
    );

    // 4. Print the result
    println!("The area under the B-spline curve (∫y dx) is: {}", area);

    // For verification, let's plot the original curve
    #[cfg(feature = "plot")]
    {
        let t = linspace(t_start, t_end, 200);
        let (x, y): (Vec<f64>, Vec<f64>) = spline.eval_vec(&t).into_iter().unzip();

        let mut plt = Plot2D::new();
        plt.set_title(&format!("Original B-Spline (Area approx. {:.4})", area))
            .set_xlabel("x")
            .set_ylabel("y")
            .insert_pair((x.clone(), y.clone()))
            .insert_pair((
                control_points.iter().map(|p| p[0]).collect(),
                control_points.iter().map(|p| p[1]).collect(),
            ))
            .set_plot_type(vec![(0, PlotType::Line), (1, PlotType::Scatter)])
            .set_color(vec![(0, "blue"), (1, "red")])
            .set_legend(vec!["Spline", "Control Points"])
            .set_style(PlotStyle::Nature)
            .set_path("example_data/bspline_with_area.png")
            .set_dpi(600)
            .savefig()?;

        println!("Generated plot: example_data/bspline_with_area.png");
    }

    Ok(())
}
