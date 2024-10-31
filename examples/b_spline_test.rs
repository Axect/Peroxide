use peroxide::fuga::*;

fn main() -> Result<(), Box<dyn Error>> {
    let knots = vec![0f64, 1f64, 2f64, 3f64];
    let degree = 3;
    let control_points = vec![
        vec![0f64, 2f64],
        vec![0.2, -1f64],
        vec![0.4, 1f64],
        vec![0.6, -1f64],
        vec![0.8, 1f64],
        vec![1f64, 2f64],
    ];

    let spline = BSpline::clamped(degree, knots, control_points.clone())?;
    let t = linspace(0f64, 3f64, 200);
    let (x, y): (Vec<f64>, Vec<f64>) = spline.eval_vec(&t).into_iter().unzip();

    #[cfg(feature = "plot")]
    {
        let control_x = control_points.iter().map(|v| v[0]).collect::<Vec<f64>>();
        let control_y = control_points.iter().map(|v| v[1]).collect::<Vec<f64>>();

        let mut plt = Plot2D::new();
        plt.insert_pair((x.clone(), y.clone()))
            .insert_pair((control_x.clone(), control_y.clone()))
            .set_plot_type(vec![(1, PlotType::Scatter)])
            .set_color(vec![(0, "darkblue"), (1, "red")])
            .set_xlabel(r"$x$")
            .set_ylabel(r"$y$")
            .set_style(PlotStyle::Nature)
            .set_dpi(600)
            .set_path("example_data/b_spline_test.png")
            .savefig()?;
    }

    let mut df = DataFrame::new(vec![]);
    df.push("t", Series::new(t));
    df.push("x", Series::new(x));
    df.push("y", Series::new(y));
    df.print();

    Ok(())
}
