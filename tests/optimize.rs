extern crate peroxide;
use peroxide::{fuga::*, hstack};

#[test]
fn test_LM() {
    let x = seq(0, 10, 0.1);
    let p_true = vec![1.0, 2.0, 3.0];
    let y = x.fmap(|t| p_true[0] * t.powi(2) + p_true[1] * t + p_true[2]);

    let p_init = vec![1f64, 1f64, 1f64];
    let data = hstack!(x, y);
    let mut opt = Optimizer::new(data, f);
    let p_est = opt
        .set_init_param(p_init)
        .set_max_iter(50)
        .set_method(LevenbergMarquardt)
        .set_lambda_init(1e-3)
        .set_lambda_max(1e+3)
        .optimize();

    p_est.print();
}

#[test]
fn test_GD() {
    let x = seq(0, 10, 0.1);
    let p_true = vec![1.0, 2.0, 3.0];
    let y = x.fmap(|t| p_true[0] * t.powi(2) + p_true[1] * t + p_true[2]);

    let p_init = vec![1f64, 1f64, 1f64];
    let data = hstack!(x, y);
    let mut opt = Optimizer::new(data, f);
    let p_est = opt
        .set_init_param(p_init)
        .set_max_iter(1000)
        .set_method(GradientDescent)
        .set_lr(1e-6)
        .optimize();

    p_est.print();
}

fn f(x: &Vec<f64>, p: Vec<AD>) -> Option<Vec<AD>> {
    Some (
        x.iter()
            .map(|t| AD1(*t, 0f64))
            .map(|t| p[0] * t.powi(2) + p[1] * t + p[2])
            .collect()
    )
}