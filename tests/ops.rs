extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_combination() {
    assert_eq!(C(10, 0), 1);
    assert_eq!(C(10, 9), 10);
    assert_eq!(C(10, 10), 1);
    assert_eq!(C(52, 5), 2598960);
    // overflowed u64 through P(30, 15) before (returned 14052247 in release)
    assert_eq!(C(30, 15), 155117520);
    assert_eq!(C(61, 30), 232714176627630544);
}

#[test]
fn test_combination_r_greater_than_n() {
    // C(n, n - r) with r > n underflowed usize -> infinite recursion in release
    assert_eq!(C(10, 12), 0);
    assert_eq!(C(0, 1), 0);
    assert_eq!(C(5, 100), 0);
}

#[test]
fn test_permutation() {
    assert_eq!(P(5, 0), 1);
    assert_eq!(P(5, 3), 60);
    assert_eq!(P(5, 5), 120);
    assert_eq!(P(30, 10), 109027350432000);
    assert_eq!(P(20, 20), 2432902008176640000);
    assert_eq!(P(3, 5), 0);
}
