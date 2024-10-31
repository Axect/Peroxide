extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_map() {
    let a = Series::new(vec![1, 2, 3, 4]);
    let b = a.map(|x: i32| x + 1);

    assert_eq!(b, Series::new(vec![2, 3, 4, 5]));
}

#[test]
fn test_mut_map() {
    let mut a = Series::new(vec![1, 2, 3, 4]);
    a.mut_map(|x: &mut i32| *x += 1);

    assert_eq!(a, Series::new(vec![2, 3, 4, 5]));
}
