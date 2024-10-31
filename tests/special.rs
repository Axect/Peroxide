use peroxide::fuga::{LambertWAccuracyMode::*, *};

#[test]
fn lambert_w_test() {
    assert_eq!(lambert_w0(1.0, Precise), 0.567143290409784);
    assert!(nearly_eq(lambert_w0(1.0, Simple), 0.567143290409784));
}
