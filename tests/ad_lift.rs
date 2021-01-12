extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_lift_ad1() {
    let a1 = AD1(2f64, 1f64);
    let x = f_ad(a1);
    x.print();
    assert_eq!(x, a1 * 2f64);
}

#[test]
fn test_lift_ad2() {
    let a2 = AD2(2f64, 1f64, 0f64);
    let x = f_ad(a2);
    x.print();
    assert_eq!(x, a2 * 2f64);
}

#[test]
fn test_lift_f64() {
    let f = 2f64;
    let lift = ADLift::new(f_ad);
    let x = lift.f_0(f);
    x.print();
    assert_eq!(x, f * 2f64);
}

fn f_ad(a: AD) -> AD {
    a * 2f64
}
