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

    let spline = BSpline::clamped(degree, knots.clone(), control_points.clone())?;

    // 2. Get the derivative and integral
    let deriv_spline = spline.derivative()?;
    let integ_spline = spline.integral()?;

    // 3. Evaluate points on all splines
    let t = linspace(0f64, 4f64, 200);
    let (x, y): (Vec<f64>, Vec<f64>) = spline.eval_vec(&t).into_iter().unzip();
    let (dx, dy): (Vec<f64>, Vec<f64>) = deriv_spline.eval_vec(&t).into_iter().unzip();
    let (ix, iy): (Vec<f64>, Vec<f64>) = integ_spline.eval_vec(&t).into_iter().unzip();

    // 4. Plotting
    #[cfg(feature = "plot")]
    {
        // Plot Original Spline and its control points
        let control_x: Vec<f64> = control_points.iter().map(|p| p[0]).collect();
        let control_y: Vec<f64> = control_points.iter().map(|p| p[1]).collect();

        let mut plt = Plot2D::new();
        plt.set_title("Original B-Spline")
            .set_xlabel("x")
            .set_ylabel("y")
            .insert_pair((x.clone(), y.clone()))
            .insert_pair((control_x, control_y))
            .set_plot_type(vec![(0, PlotType::Line), (1, PlotType::Scatter)])
            .set_color(vec![(0, "darkblue"), (1, "darkred")])
            .set_legend(vec!["Spline", "Control Points"])
            .set_style(PlotStyle::Nature)
            .set_path("example_data/bspline_original.png")
            .set_dpi(600)
            .savefig()?;

        // Plot Derivative (Hodograph) and its control points
        let deriv_control_x: Vec<f64> = deriv_spline.control_points.iter().map(|p| p[0]).collect();
        let deriv_control_y: Vec<f64> = deriv_spline.control_points.iter().map(|p| p[1]).collect();

        let mut plt_deriv = Plot2D::new();
        plt_deriv
            .set_title("B-Spline Derivative (Hodograph)")
            .set_xlabel("dx/dt")
            .set_ylabel("dy/dt")
            .insert_pair((dx, dy))
            .insert_pair((deriv_control_x, deriv_control_y))
            .set_plot_type(vec![(0, PlotType::Line), (1, PlotType::Scatter)])
            .set_color(vec![(0, "darkblue"), (1, "darkred")])
            .set_legend(vec!["Derivative", "Derivative Control Points"])
            .set_style(PlotStyle::Nature)
            .set_path("example_data/bspline_derivative.png")
            .set_dpi(600)
            .savefig()?;

        // Plot Integral and its control points
        let integ_control_x: Vec<f64> = integ_spline.control_points.iter().map(|p| p[0]).collect();
        let integ_control_y: Vec<f64> = integ_spline.control_points.iter().map(|p| p[1]).collect();

        let mut plt_integ = Plot2D::new();
        plt_integ
            .set_title("B-Spline Integral")
            .set_xlabel("Integral of x")
            .set_ylabel("Integral of y")
            .insert_pair((ix, iy))
            .insert_pair((integ_control_x, integ_control_y))
            .set_plot_type(vec![(0, PlotType::Line), (1, PlotType::Scatter)])
            .set_color(vec![(0, "darkblue"), (1, "darkred")])
            .set_legend(vec!["Integral", "Integral Control Points"])
            .set_style(PlotStyle::Nature)
            .set_path("example_data/bspline_integral.png")
            .set_dpi(600)
            .savefig()?;
    }

    println!("Generated plots for original, derivative, and integral splines.");
    println!("Please check 'example_data/bspline_original.png', 'example_data/bspline_derivative.png', and 'example_data/bspline_integral.png'.");

    Ok(())
}
