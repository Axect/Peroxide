use peroxide::fuga::{LambertWAccuracyMode::*, *};
use std::f64::consts::PI;

#[test]
fn lambert_w_test() {
    assert_eq!(lambert_w0(1.0, Precise), 0.567143290409784);
    assert!(nearly_eq(lambert_w0(1.0, Simple), 0.567143290409784));
}

#[test]
fn test_gamma_poles_and_undefined() {
    // Gamma(0) approaches infinity
    assert!(gamma(0.0).is_infinite());
    assert!(gamma(0.0).is_sign_positive());

    // Gamma for negative integers is mathematically undefined (diverges)
    assert!(gamma(-1.0).is_nan());
    assert!(gamma(-2.0).is_nan());
    assert!(gamma(-10.0).is_nan());

    // Log-Gamma goes to positive infinity for all poles
    assert!(ln_gamma(0.0).is_infinite());
    assert!(ln_gamma(-1.0).is_infinite());
    assert!(ln_gamma(-10.0).is_infinite());
}

#[test]
fn test_gamma_integer_fast_path() {
    // Standard small factorials: Gamma(n) = (n-1)!
    assert_eq!(gamma(1.0), 1.0); // 0!
    assert_eq!(gamma(2.0), 1.0); // 1!
    assert_eq!(gamma(4.0), 6.0); // 3!
    assert_eq!(gamma(5.0), 24.0); // 4!
    assert_eq!(gamma(10.0), 362_880.0); // 9!

    // Wolfram Alpha high-precision check (21!)
    // f64 can exactly represent this without precision loss
    assert_eq!(gamma(22.0), 51_090_942_171_709_440_000.0);

    // Maximum limit of f64 float representation (~171.6)
    // Ensure it doesn't panic on overflow, but correctly yields Infinity
    assert!(gamma(172.0).is_infinite());
}

#[test]
fn test_gamma_positive_floats() {
    let sqrt_pi = PI.sqrt();

    // Gamma(0.5) = sqrt(PI)
    assert!(nearly_eq(gamma(0.5), sqrt_pi));

    // Gamma(1.5) = 0.5 * sqrt(PI)
    assert!(nearly_eq(gamma(1.5), 0.5 * sqrt_pi));

    // Gamma(2.5) = 1.329340388179...
    assert!(nearly_eq(gamma(2.5), 0.75 * sqrt_pi));
}

#[test]
fn test_gamma_negative_floats_reflection() {
    let sqrt_pi = PI.sqrt();

    // Gamma(-0.5) = -2 * sqrt(PI)
    // This validates that .abs() is NOT used on the sine in gamma_approx
    assert!(nearly_eq(gamma(-0.5), -2.0 * sqrt_pi));
    assert!(gamma(-0.5).is_sign_negative());

    // Gamma(-1.5) = (4/3) * sqrt(PI)
    assert!(nearly_eq(gamma(-1.5), (4.0 / 3.0) * sqrt_pi));
    assert!(gamma(-1.5).is_sign_positive());

    // Gamma(-2.5) = -(8/15) * sqrt(PI)
    assert!(nearly_eq(gamma(-2.5), -(8.0 / 15.0) * sqrt_pi));
    assert!(gamma(-2.5).is_sign_negative());
}

#[test]
fn test_ln_gamma_consistency() {
    // ln_gamma(x) should equal ln(|Gamma(x)|) across the board
    let test_values = vec![0.5, 1.5, 2.5, 10.5];

    for &val in &test_values {
        let expected = gamma(val).ln();
        let actual = ln_gamma(val);
        assert!(
            nearly_eq(expected, actual),
            "Failed at positive float: val={}, expected={}, actual={}",
            val,
            expected,
            actual
        );
    }

    // Test Negative Floats to ensure `.abs()` prevents NaN
    let negative_test_values = vec![-0.5, -1.5, -2.5, -10.5];

    for &val in &negative_test_values {
        let expected = gamma(val).abs().ln();
        let actual = ln_gamma(val);
        assert!(
            nearly_eq(expected, actual),
            "Failed at negative float: val={}, expected={}, actual={}",
            val,
            expected,
            actual
        );
    }
}
