#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn series_map_test() {
    // Int
    let a: Vec<i32> = vec![1,2,3,4,5];
    let b: Vec<i32> = a.iter().map(|x| 2 * *x + 1).collect();
    let s_a = Series::new(a);
    let s_b = Series::new(b);
    assert_eq!(s_b, s_a.map(|t: i32| 2 * t + 1));

    // F64
    let a: Vec<f64> = vec![1.0,2.0,3.0,4.0,5.0];
    let b: Vec<f64> = a.iter().map(|x| 2f64 * *x + 1f64).collect();
    let s_a = Series::new(a);
    let s_b = Series::new(b);
    assert_eq!(s_b, s_a.map(|t: f64| 2f64 * t + 1f64));
}

#[test]
fn series_fold_test() {
    // Int
    let a: Vec<i32> = vec![1,2,3,4,5];
    let b = a.iter().fold(0, |x, y| x + *y);
    let s_a = Series::new(a);
    assert_eq!(b, s_a.fold(0i32, |x, y| x + y));

    // F64
    let a: Vec<f64> = c!(1,2,3,4,5);
    let b = a.iter().fold(0f64, |x, y| x + *y);
    let s_a = Series::new(a);
    assert_eq!(b, s_a.fold(0f64, |x, y| x + y));
}

#[test]
fn series_filter_test() {
    // Int
    let a: Vec<i32> = vec![1,2,3,4,5];
    let b: Vec<i32> = a.clone().into_iter().filter(|x| *x % 2 == 0).collect();
    let s_a = Series::new(a);
    let s_b = Series::new(b);
    assert_eq!(s_b, s_a.filter(|x: &i32| *x % 2 == 0));

    // F64
    let a: Vec<f64> = c!(1,2,3,4,5);
    let b: Vec<f64> = a.clone().into_iter().filter(|x| *x % 2f64 == 0f64).collect();
    let s_a = Series::new(a);
    let s_b = Series::new(b);
    assert_eq!(s_b, s_a.filter(|x: &f64| *x % 2f64 == 0f64));
}

#[test]
fn series_zip_with_test() {
    // Int
    let a: Vec<i32> = vec![1,2,3,4,5];
    let b: Vec<i32> = vec![5,4,3,2,1];
    let c: Vec<i32> = a.iter().zip(b.iter()).map(|(x,y)| x + y).collect();
    let s_a = Series::new(a);
    let s_b = Series::new(b);
    let s_c = Series::new(c);
    assert_eq!(s_c, s_a.zip_with(|x: i32, y: i32| x + y, &s_b));
}
