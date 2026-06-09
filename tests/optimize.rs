extern crate peroxide;
use peroxide::{fuga::*, hstack};

const P_TRUE: [f64; 3] = [1.0, 2.0, 3.0];

// y = p[0] * t^2 + p[1] * t + p[2]
#[allow(clippy::ptr_arg)]
fn quadratic(x: &Vec<f64>, p: Vec<AD>) -> Option<Vec<AD>> {
    Some(
        x.iter()
            .map(|t| AD1(*t, 0f64))
            .map(|t| p[0] * t.powi(2) + p[1] * t + p[2])
            .collect(),
    )
}

// y = p[0] * exp(p[1] * t)
#[allow(clippy::ptr_arg)]
fn exp_decay(x: &Vec<f64>, p: Vec<AD>) -> Option<Vec<AD>> {
    Some(
        x.iter()
            .map(|t| AD1(*t, 0f64))
            .map(|t| p[0] * (p[1] * t).exp())
            .collect(),
    )
}

#[test]
#[allow(non_snake_case)]
fn test_LM_quadratic_recovers_parameters() {
    let x = seq(0, 10, 0.1);
    let y = x.fmap(|t| P_TRUE[0] * t.powi(2) + P_TRUE[1] * t + P_TRUE[2]);
    let data = hstack!(x, y);

    let mut opt = Optimizer::new(data, quadratic);
    let p = opt
        .set_init_param(vec![1f64, 1f64, 1f64])
        .set_max_iter(50)
        .set_method(LevenbergMarquardt)
        .set_lambda_init(1e-3)
        .set_lambda_max(1e+3)
        .optimize();

    for i in 0..3 {
        assert!((p[i] - P_TRUE[i]).abs() < 1e-9);
    }
}

#[test]
#[allow(non_snake_case)]
fn test_GD_quadratic_makes_progress() {
    let x = seq(0, 10, 0.1);
    let y = x.fmap(|t| P_TRUE[0] * t.powi(2) + P_TRUE[1] * t + P_TRUE[2]);
    let data = hstack!(x, y);
    let p_init = [1f64, 1f64, 1f64];

    let mut opt = Optimizer::new(data, quadratic);
    let p = opt
        .set_init_param(p_init.to_vec())
        .set_max_iter(1000)
        .set_method(GradientDescent)
        .set_lr(1e-6)
        .optimize();

    let init_total: f64 = (0..3).map(|i| (p_init[i] - P_TRUE[i]).abs()).sum();
    let final_total: f64 = (0..3).map(|i| (p[i] - P_TRUE[i]).abs()).sum();
    assert!(final_total < init_total);
    for i in [1usize, 2] {
        assert!((p[i] - P_TRUE[i]).abs() < (p_init[i] - P_TRUE[i]).abs());
    }
}

#[test]
#[allow(non_snake_case)]
fn test_LM_exponential_recovers_parameters() {
    let amp = 2.0_f64;
    let rate = -0.5_f64;
    let x = seq(0, 5, 0.05);
    let y = x.fmap(|t| amp * (rate * t).exp());
    let data = hstack!(x, y);

    let mut opt = Optimizer::new(data, exp_decay);
    let p = opt
        .set_init_param(vec![1.0_f64, -0.1_f64])
        .set_max_iter(100)
        .set_method(LevenbergMarquardt)
        .set_lambda_init(1e-3)
        .set_lambda_max(1e+3)
        .optimize();

    assert!((p[0] - amp).abs() < 1e-6);
    assert!((p[1] - rate).abs() < 1e-6);
}
