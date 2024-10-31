extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_type_cast() {
    let mut a = DataFrame::new(vec![]);
    a.push("x", Series::new(vec![1, 2, 3, 4]));
    a.push("y", Series::new(vec![true, false, false, true]));

    let mut b = DataFrame::new(vec![]);
    b.push("x", Series::new(vec![1usize, 2, 3, 4]));
    b.push("y", Series::new(vec![true, false, false, true]));

    a.as_types(vec![USIZE, U8]);
    a.as_types(vec![USIZE, Bool]);

    assert_eq!(a, b);
}
