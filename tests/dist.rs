extern crate peroxide;
use peroxide::fuga::*;

#[test]
fn test_binomial() {
    let b = Binomial(100, 0.8);
    b.sample(10).print();
    assert!(nearly_eq(b.mean(), 80f64));
    assert!(nearly_eq(b.var(), 16f64));
}

#[test]
fn test_dirichlet() {
    let dir = MVDist::Dirichlet(vec![1.0, 2.0, 3.0]);

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
