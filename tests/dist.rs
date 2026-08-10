#![cfg(feature = "rand")]

extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_binomial() {
    let b = Binomial(100, 0.8);
    b.sample(10).print();
    assert!(nearly_eq(b.mean(), 80f64));
    assert!(nearly_eq(b.var(), 16f64));
}

// Reference values below cross-checked against scipy 1.17.1

#[test]
fn test_bernoulli_pdf_support() {
    let b = Bernoulli(0.3);
    assert!(nearly_eq(b.pdf(1.0), 0.3));
    assert!(nearly_eq(b.pdf(0.0), 0.7));
    // off-support values returned 1 - p before
    assert_eq!(b.pdf(-1.0), 0f64);
    assert_eq!(b.pdf(0.5), 0f64);
    assert_eq!(b.pdf(2.0), 0f64);
}

#[test]
fn test_studentt_cdf_negative_half() {
    let t3 = StudentT(3f64);
    // negative branch used cdf(-x) - 0.5 instead of 1 - cdf(-x)
    assert!(nearly_eq(t3.cdf(-5.0), 0.007696219036651147));
    assert!(nearly_eq(t3.cdf(-2.0), 0.06966298427942152));
    assert!(nearly_eq(t3.cdf(-1.0), 0.19550110947788524));
    assert!(nearly_eq(t3.cdf(-0.001), 0.49963244748473046));
    assert!(nearly_eq(t3.cdf(0.0), 0.5));
    assert!(nearly_eq(t3.cdf(1.0), 1.0 - t3.cdf(-1.0)));
    // StudentT(1) is Cauchy: cdf(-1) = 1/4 exactly
    assert!(nearly_eq(StudentT(1f64).cdf(-1.0), 0.25));
}

#[test]
fn test_binomial_pdf() {
    let b = Binomial(10, 0.3);
    assert!(nearly_eq(b.pdf(0.0), 0.0282475249));
    assert!(nearly_eq(b.pdf(3.0), 0.2668279320));
    assert!(nearly_eq(b.pdf(10.0), 5.9049e-6));
    // pdf(x) for x > n looped forever (C(n, n - r) usize underflow)
    assert_eq!(b.pdf(12.0), 0f64);
    assert_eq!(b.pdf(-1.0), 0f64);
    assert_eq!(b.pdf(5.5), 0f64);
    // C(30, 15) overflowed u64: pdf was 0.0131 in release, panic in debug
    assert!(nearly_eq(Binomial(30, 0.5).pdf(15.0), 0.14446444809436781));
    // far past u64 range
    assert!(nearly_eq(Binomial(100, 0.5).pdf(50.0), 0.07958923738717871));
    // degenerate p
    assert!(nearly_eq(Binomial(5, 0.0).pdf(0.0), 1.0));
    assert!(nearly_eq(Binomial(5, 1.0).pdf(5.0), 1.0));
}

#[test]
fn test_binomial_cdf_domain() {
    let b = Binomial(10, 0.3);
    assert!(nearly_eq(b.cdf(3.0), 0.6496107184));
    // k = n, k > n, k < 0 panicked in inc_beta before
    assert!(nearly_eq(b.cdf(10.0), 1.0));
    assert!(nearly_eq(b.cdf(11.0), 1.0));
    assert_eq!(b.cdf(-1.0), 0f64);
    // fractional k floors (was interpolated)
    assert!(nearly_eq(b.cdf(5.5), 0.9526510126));
    assert!(nearly_eq(Binomial(100, 0.5).cdf(50.0), 0.5397946186935891));
}

#[test]
fn test_dirichlet() {
    let dir = MVDist::Dirichlet(vec![1.0, 2.0, 3.0]);
    dir.sample(10).print();

    let m = dir.mean();
    assert!(nearly_eq(m[0], 1.0 / 6.0));
    assert!(nearly_eq(m[1], 1.0 / 3.0));
    assert!(nearly_eq(m[2], 0.5));

    let v = dir.var();
    assert!(nearly_eq(v[0], 5.0 / 252.0)); // 1 * 5 / (36 * 7)
    assert!(nearly_eq(v[1], 8.0 / 252.0)); // 2 * 4 / (36 * 7)
    assert!(nearly_eq(v[2], 9.0 / 252.0)); // 3 * 3 / (36 * 7)

    let pdf_val = dir.pdf(&[0.33333, 0.33333, 0.33333]);
    assert!(nearly_eq(pdf_val, 2.222155556222205));
}
