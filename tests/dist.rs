#![cfg(feature = "rand")]

use peroxide::c;
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

#[test]
fn test_wishart() {
    #[cfg(feature = "O3")]
    {
        // Define Scale Matrix V:
        // [ 1.0  0.5 ]
        // [ 0.5  1.0 ]
        let scale = matrix(c!(1.0, 0.5, 0.5, 1.0), 2, 2, Row);

        // Wishart with 5 degrees of freedom
        let w = MatDist::Wishart(5.0, scale);

        // Test Sampling
        let samples = w.sample(10);
        samples.print();
        assert_eq!(samples.len(), 10);
        assert_eq!(samples[0].nrow(), 2);
        assert_eq!(samples[0].ncol(), 2);

        // Test Mean (E[X] = df * V)
        let m = w.mean();
        assert!(nearly_eq(m[(0, 0)], 5.0));
        assert!(nearly_eq(m[(0, 1)], 2.5));
        assert!(nearly_eq(m[(1, 0)], 2.5));
        assert!(nearly_eq(m[(1, 1)], 5.0));

        // Test Variance (Element-wise: Var(X_ij) = df * (V_ij^2 + V_ii * V_jj))
        let v = w.var();
        assert!(nearly_eq(v[(0, 0)], 10.0)); // 5 * (1^2 + 1*1)
        assert!(nearly_eq(v[(0, 1)], 6.25)); // 5 * (0.5^2 + 1*1) = 5 * 1.25
        assert!(nearly_eq(v[(1, 0)], 6.25));
        assert!(nearly_eq(v[(1, 1)], 10.0));

        // Test Covariance (Vectorized 4x4 matrix)
        // Cov(X_ij, X_kl) = df * (V_ik * V_jl + V_il * V_jk)
        let cov = w.cov();
        assert_eq!(cov.nrow(), 4);
        assert_eq!(cov.ncol(), 4);

        // Cov(X_00, X_00) -> mapped to index (0, 0)
        assert!(nearly_eq(cov[(0, 0)], 10.0));

        // Cov(X_00, X_11) -> mapped to index (0, 3)
        // 5 * (V_01 * V_01 + V_01 * V_01) = 5 * (0.25 + 0.25) = 2.5
        assert!(nearly_eq(cov[(0, 3)], 2.5));

        // Cov(X_01, X_01) -> mapped to index (1, 1)
        // 5 * (V_00 * V_11 + V_01 * V_10) = 5 * (1 + 0.25) = 6.25
        assert!(nearly_eq(cov[(1, 1)], 6.25));

        // Cov(X_01, X_10) -> mapped to index (1, 2)
        // 5 * (V_01 * V_10 + V_00 * V_11) = 5 * (0.25 + 1.0) = 6.25
        assert!(nearly_eq(cov[(1, 2)], 6.25));

        // 5. Test PDF against an exact closed-form calculation
        // Let test_x = [ 2.0  0.0 ]
        //              [ 0.0  2.0 ]
        let test_x = matrix(c!(2.0, 0.0, 0.0, 2.0), 2, 2, Row);
        let pdf_val = w.pdf(&test_x);

        // Exact closed-form math for this specific test case:
        // Numerator:   |X|^1 * exp(-0.5 * tr(V^-1 * X)) = 4 * exp(-8/3)
        // Denominator: 2^5 * |V|^2.5 * Gamma_2(2.5) = 32 * (0.75)^2.5 * (0.75 * pi)
        // Which simplifies cleanly to: 16 * exp(-8/3) / (27 * sqrt(3) * pi)
        let expected_pdf =
            16.0 * (-8.0f64 / 3.0).exp() / (27.0 * 3.0f64.sqrt() * std::f64::consts::PI);

        assert!(nearly_eq(pdf_val, expected_pdf));
    }
}

#[test]
fn test_wishart_1d_chi_square() {
    #[cfg(feature = "O3")]
    {
        let df = 5.0;
        let scale = matrix(c!(1.0), 1, 1, Row);
        let w = MatDist::Wishart(df, scale);

        // A Chi-Square(df) distribution has Mean = df, Variance = 2*df
        let m = w.mean();
        let v = w.var();

        assert!(nearly_eq(m[(0, 0)], df));
        assert!(nearly_eq(v[(0, 0)], 2.0 * df));

        // Check PDF Convergence
        let x_val = 3.0;
        let test_x = matrix(c!(x_val), 1, 1, Row);
        let wishart_pdf = w.pdf(&test_x);

        let chi_pdf = (x_val.powf(df / 2.0 - 1.0) * (-x_val / 2.0).exp())
            / (2.0f64.powf(df / 2.0) * gamma(df / 2.0));

        assert!(nearly_eq(wishart_pdf, chi_pdf));
    }
}

#[test]
fn test_wishart_support_and_invalid_inputs() {
    let scale = matrix(c!(1.0, 0.5, 0.5, 1.0), 2, 2, Row);
    let w = MatDist::Wishart(5.0, scale);

    // x = -I (Negative Definite) -> ln_pdf should be -inf, pdf should be 0
    let test_x_neg = matrix(c!(-1.0, 0.0, 0.0, -1.0), 2, 2, Row);
    assert_eq!(w.ln_pdf(&test_x_neg), std::f64::NEG_INFINITY);
    assert_eq!(w.pdf(&test_x_neg), 0.0);

    // x = Asymmetric -> ln_pdf should be -inf
    let test_x_asym = matrix(c!(1.0, 0.8, 0.2, 1.0), 2, 2, Row);
    assert_eq!(w.ln_pdf(&test_x_asym), std::f64::NEG_INFINITY);
}

#[test]
#[should_panic]
fn test_wishart_invalid_df_panic() {
    // df = 0.5, but p = 2. This is invalid and should panic
    let scale = matrix(c!(1.0, 0.0, 0.0, 1.0), 2, 2, Row);
    let w = MatDist::Wishart(0.5, scale);

    // We only need a valid SPD matrix to pass the bounds check and trigger the df panic
    let test_x = matrix(c!(1.0, 0.0, 0.0, 1.0), 2, 2, Row);
    w.ln_pdf(&test_x);
}

#[test]
fn test_wishart_reproducibility() {
    #[cfg(feature = "O3")]
    {
        let mut rng1 = smallrng_from_seed(42);
        let mut rng2 = smallrng_from_seed(42);

        let scale = matrix(c!(1.0, 0.5, 0.5, 1.0), 2, 2, Row);
        let w = MatDist::Wishart(5.0, scale);

        let s1 = w.sample_with_rng(&mut rng1, 1);
        let s2 = w.sample_with_rng(&mut rng2, 1);

        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(s1[0][(i, j)], s2[0][(i, j)]);
            }
        }
    }
}

#[test]
fn test_wishart_empirical_convergence() {
    #[cfg(feature = "O3")]
    {
        let mut rng = smallrng_from_seed(42);
        let scale = matrix(c!(1.0, 0.5, 0.5, 1.0), 2, 2, Row);
        let w = MatDist::Wishart(5.0, scale);

        let n_samples = 20_000;
        let samples = w.sample_with_rng(&mut rng, n_samples);

        let p = 2;
        let p_sq = p * p;
        let n_f64 = n_samples as f64;

        let mut emp_mean = zeros(p, p);
        for s in &samples {
            emp_mean = &emp_mean + s;
        }

        for i in 0..p {
            for j in 0..p {
                emp_mean[(i, j)] /= n_f64;
            }
        }

        let mut emp_cov = zeros(p_sq, p_sq);
        for s in &samples {
            for i in 0..p {
                for j in 0..p {
                    for k in 0..p {
                        for l in 0..p {
                            let row_idx = i * p + j;
                            let col_idx = k * p + l;

                            let diff1 = s[(i, j)] - emp_mean[(i, j)];
                            let diff2 = s[(k, l)] - emp_mean[(k, l)];

                            emp_cov[(row_idx, col_idx)] += diff1 * diff2;
                        }
                    }
                }
            }
        }

        for i in 0..p_sq {
            for j in 0..p_sq {
                emp_cov[(i, j)] /= n_f64 - 1.0;
            }
        }

        let theo_mean = w.mean();
        let theo_cov = w.cov();

        // Ensure empirical mean is within a small epsilon of theoretical mean
        for i in 0..p {
            for j in 0..p {
                assert!(
                    (emp_mean[(i, j)] - theo_mean[(i, j)]).abs() < 0.1,
                    "Mean failed to converge at ({}, {})",
                    i,
                    j
                );
            }
        }

        // We have a slightly higher epsilon as covariance uses higher moments
        for i in 0..p_sq {
            for j in 0..p_sq {
                assert!(
                    (emp_cov[(i, j)] - theo_cov[(i, j)]).abs() < 0.5,
                    "Covariance failed to converge at ({}, {})",
                    i,
                    j
                );
            }
        }
    }
}
