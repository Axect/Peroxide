use std::f64::consts::{E, FRAC_PI_2, PI};

use peroxide::fuga::*;

const TOL: f64 = 1e-10;

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < TOL,
        "actual = {:.17e}, expected = {:.17e}, diff = {:.3e}",
        actual,
        expected,
        (actual - expected).abs()
    );
}

fn gk_methods(tol: f64) -> [Integral; 6] {
    [
        G7K15R(tol, 20),
        G10K21R(tol, 20),
        G15K31R(tol, 20),
        G20K41R(tol, 20),
        G25K51R(tol, 20),
        G30K61R(tol, 20),
    ]
}

#[test]
fn test_polynomial() {
    // int_0^1 x^3 dx = 1/4
    for m in gk_methods(1e-8) {
        let r: f64 = integrate(|x: f64| x.powi(3), (0.0, 1.0), m);
        assert_close(r, 0.25);
    }
}

#[test]
fn test_trig() {
    // int_0^1 cos x dx = sin(1)
    for m in gk_methods(1e-8) {
        let r: f64 = integrate(|x: f64| x.cos(), (0.0, 1.0), m);
        assert_close(r, 1f64.sin());
    }
    // int_0^pi sin x dx = 2
    for m in gk_methods(1e-8) {
        let r: f64 = integrate(|x: f64| x.sin(), (0.0, PI), m);
        assert_close(r, 2.0);
    }
}

#[test]
fn test_exponential() {
    // int_0^2 exp(x) dx = e^2 - 1
    for m in gk_methods(1e-8) {
        let r: f64 = integrate(|x: f64| x.exp(), (0.0, 2.0), m);
        assert_close(r, E * E - 1.0);
    }
}

#[test]
fn test_gaussian() {
    // int_0^1 exp(-x^2) dx = (sqrt(pi)/2) * erf(1) = 0.746824132812427...
    let expected = 0.746_824_132_812_427_f64;
    for m in gk_methods(1e-10) {
        let r: f64 = integrate(|x: f64| (-x * x).exp(), (0.0, 1.0), m);
        assert_close(r, expected);
    }
}

#[test]
fn test_rational_symmetric() {
    // int_{-1}^{1} 1 / (1 + x^2) dx = pi/2
    for m in gk_methods(1e-10) {
        let r: f64 = integrate(|x: f64| 1.0 / (1.0 + x * x), (-1.0, 1.0), m);
        assert_close(r, FRAC_PI_2);
    }
}

#[test]
fn test_oscillatory() {
    // int_0^{2 pi} sin(5 x) dx = 0
    // Target value is zero, so the relative-tolerance variants would not
    // converge; use the absolute-tolerance ones.
    let abs_methods = [
        G7K15(1e-10, 20),
        G10K21(1e-10, 20),
        G15K31(1e-10, 20),
        G20K41(1e-10, 20),
        G25K51(1e-10, 20),
        G30K61(1e-10, 20),
    ];
    for m in abs_methods {
        let r: f64 = integrate(|x: f64| (5.0 * x).sin(), (0.0, 2.0 * PI), m);
        assert_close(r, 0.0);
    }
}

#[test]
fn test_endpoint_independence() {
    let m = || G10K21R(1e-10, 20);
    let f = |x: f64| (-x).exp() + x.sin();
    let full: f64 = integrate(f, (0.0, 3.0), m());
    let left: f64 = integrate(f, (0.0, 1.7), m());
    let right: f64 = integrate(f, (1.7, 3.0), m());
    assert_close(full, left + right);
}
