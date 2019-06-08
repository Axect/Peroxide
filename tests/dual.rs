extern crate peroxide;
use peroxide::*;

#[test]
fn test_dot_dual() {
    let x = vec![Dual::new(1,5), Dual::new(2,4), Dual::new(3,3), Dual::new(4,2)];
    let y = vec![Dual::new(5,1), Dual::new(4,2), Dual::new(3,3), Dual::new(2,4)];
    assert_eq!(Dual::new(30,84), x.dot(&y));
}