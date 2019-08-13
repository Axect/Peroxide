#[macro_use]
extern crate peroxide;
use peroxide::*;

fn main() {
    let normal = Normal(0f64, 0.1f64);
    let normal2 = Normal(0f64, 100f64);

    let mut x = seq(0., 99., 1f64);
    x = zip_with(|a, b| (a + b).abs(), &x, &normal.sample(x.len()));

    let mut y = x.fmap(|t| t.powi(2));
    y = zip_with(|a, b| a + b, &y, &normal2.sample(y.len()));

    let n_init = vec![1f64];
    let data = hstack!(x.clone(), y.clone());

    let mut opt = Optimizer::new(data, quad);
    let p = opt
        .set_init_param(n_init)
        .set_max_iter(50)
        .set_method(LevenbergMarquardt)
        .optimize();
    p.print();
    opt.get_error().print();

    let z = quad(&x, NumberVector::from_f64_vec(p)).to_f64_vec();

    let mut plt = Plot2D::new();
    plt.set_domain(x)
        .insert_image(y)
        .insert_image(z)
        .set_legend(vec!["Noise", "Fit"])
        .set_title("Test ($y=x^2$ with noise)")
        .set_path("example_data/lm_test.png")
        .set_marker(vec![Point, Line])
        .savefig()
        .expect("Can't draw a plot");
}

fn quad(x: &Vec<f64>, n: Vec<Number>) -> Vec<Number> {
    x.clone()
        .into_iter()
        .map(|t| Number::from_f64(t))
        .map(|t| t.powf(n[0]))
        .collect()
}
