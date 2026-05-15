//! Regression tests for `Optimizer` (`tests/optimize.rs`).
//!
//! The original version only called `.optimize()` and printed the
//! returned parameter vector, which would not have caught a silent
//! regression in the optimizer itself. These tests assert the
//! parameters that the fit should recover.

extern crate peroxide;
use peroxide::{fuga::*, hstack};

const P_TRUE: [f64; 3] = [1.0, 2.0, 3.0];

// Quadratic model y = p[0] * t^2 + p[1] * t + p[2].
//
// Signature must stay `&Vec<f64>` to satisfy `Optimizer`'s trait bound.
#[allow(clippy::ptr_arg)]
fn quadratic(x: &Vec<f64>, p: Vec<AD>) -> Option<Vec<AD>> {
    Some(
        x.iter()
            .map(|t| AD1(*t, 0f64))
            .map(|t| p[0] * t.powi(2) + p[1] * t + p[2])
            .collect(),
    )
}

// Exponential model y = p[0] * exp(p[1] * t).
//
// Signature must stay `&Vec<f64>` to satisfy `Optimizer`'s trait bound.
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
    let p_est = opt
        .set_init_param(vec![1f64, 1f64, 1f64])
        .set_max_iter(50)
        .set_method(LevenbergMarquardt)
        .set_lambda_init(1e-3)
        .set_lambda_max(1e+3)
        .optimize();

    // Levenberg-Marquardt on a linear-in-parameters least-squares
    // problem should converge to the exact solution within a handful
    // of iterations. A tight tolerance is appropriate here.
    for (i, expected) in P_TRUE.iter().enumerate() {
        let err = (p_est[i] - expected).abs();
        assert!(
            err < 1e-9,
            "LM did not recover parameter {i}: got {got}, expected {expected}, |err| = {err:.2e}",
            i = i,
            got = p_est[i],
            expected = expected,
            err = err,
        );
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
    let p_est = opt
        .set_init_param(p_init.to_vec())
        .set_max_iter(1000)
        .set_method(GradientDescent)
        .set_lr(1e-6)
        .optimize();

    // Vanilla gradient descent on this problem is intrinsically slow
    // because the design matrix [t^2, t, 1] over t in [0, 10] has a
    // poorly conditioned Hessian (parameter scales differ by orders
    // of magnitude). After 1000 iterations the fit has not converged
    // to high accuracy, but every parameter should be measurably
    // closer to the truth than the initial guess. This is what
    // catches a regression that breaks the GD direction.
    let initial_total: f64 = (0..3).map(|i| (p_init[i] - P_TRUE[i]).abs()).sum();
    let final_total: f64 = (0..3).map(|i| (p_est[i] - P_TRUE[i]).abs()).sum();
    assert!(
        final_total < initial_total,
        "GD did not reduce total parameter error: initial = {:.3e}, final = {:.3e}",
        initial_total,
        final_total,
    );
    // The two parameters that were not already exact should each have
    // moved toward their true values.
    for i in [1usize, 2] {
        let init_err = (p_init[i] - P_TRUE[i]).abs();
        let final_err = (p_est[i] - P_TRUE[i]).abs();
        assert!(
            final_err < init_err,
            "GD did not improve parameter {}: initial |err| = {:.3e}, final |err| = {:.3e}",
            i,
            init_err,
            final_err,
        );
    }
}

#[test]
#[allow(non_snake_case)]
fn test_LM_exponential_recovers_parameters() {
    // y = 2 * exp(-0.5 * t).
    let amp_true = 2.0_f64;
    let rate_true = -0.5_f64;
    let x = seq(0, 5, 0.05);
    let y = x.fmap(|t| amp_true * (rate_true * t).exp());
    let data = hstack!(x, y);

    let mut opt = Optimizer::new(data, exp_decay);
    let p_est = opt
        .set_init_param(vec![1.0_f64, -0.1_f64])
        .set_max_iter(100)
        .set_method(LevenbergMarquardt)
        .set_lambda_init(1e-3)
        .set_lambda_max(1e+3)
        .optimize();

    let amp_err = (p_est[0] - amp_true).abs();
    let rate_err = (p_est[1] - rate_true).abs();
    // LM converges to high accuracy on this nonlinear-in-parameters
    // model when the starting point is in the basin of attraction.
    assert!(
        amp_err < 1e-6,
        "LM did not recover amplitude: got {got}, expected {amp_true}, |err| = {amp_err:.2e}",
        got = p_est[0],
    );
    assert!(
        rate_err < 1e-6,
        "LM did not recover rate: got {got}, expected {rate_true}, |err| = {rate_err:.2e}",
        got = p_est[1],
    );
}
