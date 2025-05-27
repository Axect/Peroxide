use peroxide::fuga::*;

#[test]
fn test_cosine_integral() {
    let a = 0.0;
    let b = 1.0;
    let result = integrate(|x| x.cos(), (a, b), G7K15R(1e-8, 20));
    let expected = 1f64.sin();
    assert!(
        (result - expected).abs() < 1e-10,
        "Expected: {}, Got: {}",
        expected,
        result
    );
}
