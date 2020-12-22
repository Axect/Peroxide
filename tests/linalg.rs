extern crate peroxide;
#[allow(unused_imports)]
use peroxide::fuga::*;

#[cfg(feature = "nc")]
#[test]
fn test_inverse() {
    for i in 5..21 {
        let df = DataFrame::read_nc_by_header(
            &format!("test_data/rand_mat/randmat_{}.nc", i),
            vec!["m", "inv"],
        )
        .unwrap();
        let a_vec: Vec<f64> = df["m"].to_vec();
        let b_vec: Vec<f64> = df["inv"].to_vec();
        let a: Matrix = matrix(a_vec, i, i, Col);
        let b: Matrix = matrix(b_vec, i, i, Col);
        let c = a.inv();
        assert_eq!(b, c);
    }
}

#[cfg(feature = "nc")]
#[test]
fn test_matmul() {
    for i in 5..21 {
        let df = DataFrame::read_nc_by_header(
            &format!("test_data/rand_mat/randmat_{}.nc", i),
            vec!["m", "mm"],
        )
        .unwrap();
        let a_vec: Vec<f64> = df["m"].to_vec();
        let b_vec: Vec<f64> = df["mm"].to_vec();
        let a: Matrix = matrix(a_vec, i, i, Col);
        let b: Matrix = matrix(b_vec, i, i, Col);
        let c = &a * &a;
        assert_eq!(b, c);
    }
}

#[cfg(feature = "nc")]
#[test]
fn test_matvecmul() {
    for i in 5..21 {
        let df = DataFrame::read_nc_by_header(
            &format!("test_data/rand_mat/randmat_{}.nc", i),
            vec!["m", "v", "mv"],
        )
        .unwrap();
        let a_vec: Vec<f64> = df["m"].to_vec();
        let a: Matrix = matrix(a_vec, i, i, Col);
        let b: Vec<f64> = df["v"].to_vec();
        let c: Vec<f64> = df["mv"].to_vec();
        let d = &a * &b;
        assert!(eq_vec(&c, &d, 1e-8));
    }
}

#[cfg(feature = "nc")]
#[test]
fn test_det() {
    for i in 5..21 {
        let df = DataFrame::read_nc_by_header(
            &format!("test_data/rand_mat/randmat_{}.nc", i),
            vec!["m", "det"],
        )
        .unwrap();
        let a_vec: Vec<f64> = df["m"].to_vec();
        let a: Matrix = matrix(a_vec, i, i, Col);
        let b: f64 = df["det"].at(0).unwrap();
        let c = a.det();
        assert!((b - c).abs() <= 1e-10);
    }
}

#[cfg(feature = "nc")]
#[test]
fn test_pinv() {
    for i in 5..21 {
        let df = DataFrame::read_nc_by_header(
            &format!("test_data/rand_mat/randmat_{}.nc", i),
            vec!["pm", "pinv"],
        )
        .unwrap();
        let a_vec: Vec<f64> = df["pm"].to_vec();
        let b_vec: Vec<f64> = df["pinv"].to_vec();
        let a: Matrix = matrix(a_vec, i + 1, i - 1, Col);
        let b: Matrix = matrix(b_vec, i - 1, i + 1, Col);
        let c = a.pseudo_inv();
        assert_eq!(b, c);
    }
}

#[cfg(feature = "nc")]
#[test]
fn test_solve() {
    for i in 5..21 {
        let df = DataFrame::read_nc_by_header(
            &format!("test_data/rand_mat/randmat_{}.nc", i),
            vec!["m", "v", "sol"],
        )
        .unwrap();
        let a_vec: Vec<f64> = df["m"].to_vec();
        let a: Matrix = matrix(a_vec, i, i, Col);
        let b: Vec<f64> = df["v"].to_vec();
        let x: Vec<f64> = df["sol"].to_vec();
        let c = a.solve(&b, LU);
        let d = a.solve(&b, WAZ);
        x.print();
        c.print();
        d.print();
        assert!(eq_vec(&x, &c, 1e-6));
        assert!(eq_vec(&x, &d, 1e-6));
    }
}
