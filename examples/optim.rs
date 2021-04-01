#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    // To prepare noise
    let normal = Normal(0f64, 0.1f64);
    let normal2 = Normal(0f64, 100f64);

    // Noise to domain
    let mut x = seq(0., 99., 1f64);
    x = zip_with(|a, b| (a + b).abs(), &x, &normal.sample(x.len()));

    // Noise to image
    let mut y = x.fmap(|t| t.powi(2));
    y = zip_with(|a, b| a + b, &y, &normal2.sample(y.len()));

    // Initial parameter
    let n_init = vec![1f64];
    let data = hstack!(x.clone(), y.clone());

    // Optimizer setting
    let mut opt = Optimizer::new(data, quad);
    let p = opt.set_init_param(n_init)
        .set_max_iter(50)
        .set_method(LevenbergMarquardt)
        .optimize();
    p.print();                  // Optimized parameter
    opt.get_error().print();    // Optimized RMSE
}

fn quad(x: &Vec<f64>, n: Vec<AD>) -> Option<Vec<AD>> {
    Some(
        x.clone().into_iter()
            .map(|t| AD1(t, 0f64))
            .map(|t| t.pow(n[0]))
            .collect()
    )
}
