// #[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_horner_division() {
    let a = Polynomial::new(vec![1f64, -4f64, 4f64, 3f64, -8f64, 4f64]);
    let b = Polynomial::new(vec![1f64, -2f64]);

    let (c, remainder) = a.horner_division(&b);
    assert_eq!(c.coef, vec![1f64, -2f64, 0f64, 3f64, -2f64]);
    assert_eq!(remainder, 0f64);
}

#[test]
fn test_translate_x() {
    let a = Polynomial::new(vec![1f64, -4f64, 4f64, 3f64, -8f64, 4f64]);
    let b = a.translate_x(-6);

    for i in -10..10 {
        assert_eq!(a.eval(i), b.eval(i - 6));
    }
}