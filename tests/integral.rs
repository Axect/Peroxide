//! Regression tests for `integrate(...)` (`tests/integral.rs`).
//!
//! The original version only exercised `int x.cos() dx` over [0, 1]
//! with `G7K15R`. These tests cover several integrand families
//! (polynomial, trigonometric, exponential, Gaussian, rational) and
//! sweep every Gauss-Kronrod tolerance variant so that an accidental
//! change to any one quadrature path is caught.

use std::f64::consts::{E, PI};

use peroxide::fuga::*;

/// Tolerance for smooth, low-amplitude integrands. The adaptive
/// Gauss-Kronrod methods converge to machine epsilon on the test
/// integrands; the tolerance leaves a comfortable safety margin for
/// platform-specific rounding without masking a real regression.
const SMOOTH_TOL: f64 = 1e-10;

fn assert_close(name: &str, method: &str, got: f64, expected: f64, tol: f64) {
    let err = (got - expected).abs();
    assert!(
        err < tol,
        "{}: {} gave {:.17e}, expected {:.17e}, |err| = {:.2e}",
        name,
        method,
        got,
        expected,
        err,
    );
}

/// All adaptive Gauss-Kronrod variants exposed by the library.
fn gk_methods(tol: f64) -> [(&'static str, Integral); 6] {
    [
        ("G7K15R", G7K15R(tol, 20)),
        ("G10K21R", G10K21R(tol, 20)),
        ("G15K31R", G15K31R(tol, 20)),
        ("G20K41R", G20K41R(tol, 20)),
        ("G25K51R", G25K51R(tol, 20)),
        ("G30K61R", G30K61R(tol, 20)),
    ]
}

#[test]
fn test_polynomial() {
    // int_0^1 x^3 dx = 1/4
    let expected = 0.25_f64;
    for (name, method) in gk_methods(1e-8) {
        let got: f64 = integrate(|x: f64| x.powi(3), (0.0, 1.0), method);
        assert_close("int x^3 over [0,1]", name, got, expected, SMOOTH_TOL);
    }
}

#[test]
fn test_trig() {
    // int_0^1 cos(x) dx = sin(1)
    let expected = 1f64.sin();
    for (name, method) in gk_methods(1e-8) {
        let got: f64 = integrate(|x: f64| x.cos(), (0.0, 1.0), method);
        assert_close("int cos over [0,1]", name, got, expected, SMOOTH_TOL);
    }

    // int_0^pi sin(x) dx = 2
    let expected = 2.0_f64;
    for (name, method) in gk_methods(1e-8) {
        let got: f64 = integrate(|x: f64| x.sin(), (0.0, PI), method);
        assert_close("int sin over [0,pi]", name, got, expected, SMOOTH_TOL);
    }
}

#[test]
fn test_exponential() {
    // int_0^2 exp(x) dx = exp(2) - 1
    let expected = E * E - 1.0;
    for (name, method) in gk_methods(1e-8) {
        let got: f64 = integrate(|x: f64| x.exp(), (0.0, 2.0), method);
        assert_close("int exp over [0,2]", name, got, expected, SMOOTH_TOL);
    }
}

#[test]
fn test_gaussian() {
    // int_0^1 exp(-x^2) dx = (sqrt(pi)/2) * erf(1)
    //   = 0.74682413281242702540 (high-precision reference)
    let expected = 0.746_824_132_812_427_0_f64;
    for (name, method) in gk_methods(1e-10) {
        let got: f64 = integrate(|x: f64| (-x * x).exp(), (0.0, 1.0), method);
        assert_close(
            "int exp(-x^2) over [0,1]",
            name,
            got,
            expected,
            SMOOTH_TOL,
        );
    }
}

#[test]
fn test_rational_symmetric() {
    // int_{-1}^{1} 1 / (1 + x^2) dx = pi/2.
    // The integrand peaks at x = 0; this is a good stress for the
    // adaptive subdivision logic.
    let expected = std::f64::consts::FRAC_PI_2;
    for (name, method) in gk_methods(1e-10) {
        let got: f64 = integrate(|x: f64| 1.0 / (1.0 + x * x), (-1.0, 1.0), method);
        assert_close(
            "int 1/(1+x^2) over [-1,1]",
            name,
            got,
            expected,
            SMOOTH_TOL,
        );
    }
}

#[test]
fn test_oscillatory() {
    // int_0^{2 pi} sin(5 x) dx = 0
    //
    // Pure cancellation makes this an easy way to spot a quadrature
    // path that is mis-summing the contributions: any non-trivial
    // bias would survive the cancellation. Use the absolute-tolerance
    // variants here because a target of 0 makes the relative
    // tolerance (`tol * |I|`) collapse to zero.
    let abs_methods: [(&str, Integral); 6] = [
        ("G7K15", G7K15(1e-10, 20)),
        ("G10K21", G10K21(1e-10, 20)),
        ("G15K31", G15K31(1e-10, 20)),
        ("G20K41", G20K41(1e-10, 20)),
        ("G25K51", G25K51(1e-10, 20)),
        ("G30K61", G30K61(1e-10, 20)),
    ];
    let expected = 0.0_f64;
    for (name, method) in abs_methods {
        let got: f64 = integrate(|x: f64| (5.0 * x).sin(), (0.0, 2.0 * PI), method);
        assert_close(
            "int sin(5x) over [0,2pi]",
            name,
            got,
            expected,
            SMOOTH_TOL,
        );
    }
}

#[test]
fn test_endpoint_independence() {
    // Splitting the interval at an interior point must reproduce the
    // value of the integral over the full range, within tolerance.
    let method = || G10K21R(1e-10, 20);
    let f = |x: f64| (-x).exp() + x.sin();
    let full: f64 = integrate(f, (0.0, 3.0), method());
    let left: f64 = integrate(f, (0.0, 1.7), method());
    let right: f64 = integrate(f, (1.7, 3.0), method());
    assert!(
        (full - (left + right)).abs() < SMOOTH_TOL,
        "Split integral does not match the full range: full = {:.17e}, left + right = {:.17e}",
        full,
        left + right,
    );
}
