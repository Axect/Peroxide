extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_lift_ad2_ad1() {
    let a1 = AD1::new(2f64, 1f64);
    let lift = ADLift::new(f_ad2);
    let x = lift.call_stable(a1);
    x.print();
    assert_eq!(x, a1 * 2f64);
}

#[test]
fn test_lift_ad2_f64() {
    let f = 2f64;
    let lift = ADLift::new(f_ad2);
    let x = lift.call_stable(f);
    x.print();
    assert_eq!(x, f * 2f64);
}

fn f_ad2(a2: AD2) -> AD2 {
    a2 * 2f64
}
