extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_mean() {
    let a: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    assert_eq!(a.mean(), 3.0);
}

#[test]
fn test_mean_stable() {
    let a: Vec<f32> = vec![1.0; 10000000];
    let diff = 10000.0;
    let b = a.iter().map(|x| x + diff).collect::<Vec<f32>>();
    assert_eq!(a.mean(), b.mean() - diff);
}

#[test]
fn test_variance() {
    let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    assert_eq!(a.var(), 2.5);
}

#[test]
fn test_variance_stable() {
    let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let diff = 1000000000.0;
    let b = a.iter().map(|x| x + diff).collect::<Vec<f64>>();
    assert_eq!(a.var(), b.var());
}
