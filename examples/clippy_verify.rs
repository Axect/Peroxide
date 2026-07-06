//! Determinism oracle for the Phase 2 `needless_range_loop` refactor
//! (JOSS review feedback from ariostas).
//!
//! Each `cargo clippy --fix` suggestion for an index-based loop is
//! semantically equivalent to the original only if it preserves bitwise
//! output. This example exercises every public entry point that
//! transitively invokes one of the 30 affected loops, with seeded
//! random inputs, and prints a hash + a small spot-check value for
//! each result.
//!
//! Workflow:
//!
//! ```ignore
//! # snapshot baseline before touching algorithm code
//! cargo run --example clippy_verify > /tmp/clippy_verify.golden
//! # apply needless_range_loop refactors, then verify equivalence
//! cargo run --example clippy_verify > /tmp/clippy_verify.after
//! diff /tmp/clippy_verify.golden /tmp/clippy_verify.after  # must be empty
//! ```

use std::hash::{DefaultHasher, Hash, Hasher};

use peroxide::fuga::*;
use peroxide::numerical::integral::kronrod_quadrature;
use peroxide::numerical::interp::chebyshev_nodes;
use peroxide::special::lanczos::ln_gamma_approx;
use peroxide::structure::sparse::SPMatrix;
use peroxide::util::useful::{max, min};
use rand::{rngs::SmallRng, SeedableRng};

fn hash_f64s(v: &[f64]) -> String {
    let mut h = DefaultHasher::new();
    v.len().hash(&mut h);
    for x in v {
        x.to_bits().hash(&mut h);
    }
    format!("{:016x}", h.finish())
}

fn hash_matrix(m: &Matrix) -> String {
    let mut h = DefaultHasher::new();
    m.nrow().hash(&mut h);
    m.ncol().hash(&mut h);
    // the layout enum has no Hash impl; use Debug as a stable surrogate
    format!("{:?}", m.layout()).hash(&mut h);
    for x in m.as_slice() {
        x.to_bits().hash(&mut h);
    }
    format!("{:016x}", h.finish())
}

fn print_kv<T: std::fmt::Display>(label: &str, value: T) {
    println!("{:48} {}", label, value);
}

fn check_integral() {
    // kronrod_quadrature drives `compute_gauss_kronrod_sums_stored`,
    // which is the function holding the flagged loop. Call it directly
    // at every supported order so the index-based body fires for each
    // node count.
    let f = |x: f64| x.cos() + 0.25 * x * x;
    for &n in &[15_usize, 21, 31, 41, 51, 61] {
        let r: f64 = kronrod_quadrature(f, n, (0.0, 1.7));
        print_kv(&format!("kronrod_quadrature n={}", n), format!("{r:.17e}"));
        // also exercise with a different smooth integrand to widen coverage
        let r2: f64 = kronrod_quadrature(|x: f64| (-x * x).exp(), n, (-2.5, 1.5));
        print_kv(
            &format!("kronrod_quadrature gauss n={}", n),
            format!("{r2:.17e}"),
        );
    }
}

fn check_chebyshev_lagrange() {
    // chebyshev_nodes carries one of the affected loops directly.
    for &n in &[4_usize, 8, 16, 32] {
        let nodes = chebyshev_nodes(n, -2.0, 3.5);
        print_kv(&format!("chebyshev_nodes(n={n}) hash"), hash_f64s(&nodes));
    }
    // lagrange_polynomial uses divided-difference loops touched by clippy.
    let xs = vec![0.0_f64, 0.5, 1.1, 1.7, 2.4, 3.0];
    let ys: Vec<f64> = xs.iter().map(|x| x.sin() + 0.5 * x).collect();
    let p = lagrange_polynomial(xs.clone(), ys.clone());
    for &q in &[0.1_f64, 0.7, 1.3, 1.9, 2.5] {
        print_kv(
            &format!("lagrange.eval({q})"),
            format!("{:.17e}", p.eval(q)),
        );
    }
    print_kv(
        "lagrange.derivative coeffs",
        hash_f64s(&p.derivative().coef),
    );
    print_kv("lagrange.integral coeffs", hash_f64s(&p.integral().coef));
}

fn check_lanczos() {
    for &z in &[0.5_f64, 1.0, 1.5, 2.5, 3.7, 5.5, 10.25, 100.0] {
        print_kv(
            &format!("ln_gamma_approx({z})"),
            format!("{:.17e}", ln_gamma_approx(z)),
        );
    }
}

fn check_ode() {
    struct StiffSpiral;
    impl ODEProblem for StiffSpiral {
        fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) -> anyhow::Result<()> {
            dy[0] = -0.5 * y[0] - 1.5 * y[1];
            dy[1] = 1.5 * y[0] - 0.5 * y[1];
            Ok(())
        }
    }
    macro_rules! run_explicit {
        ($name:expr, $integ:expr) => {{
            let (t, y) = BasicODESolver::new($integ)
                .solve(&StiffSpiral, (0.0, 1.5), 0.05, &[1.0, 0.0])
                .unwrap();
            let y_flat: Vec<f64> = y.iter().flatten().copied().collect();
            print_kv(&format!("ode {} t-hash", $name), hash_f64s(&t));
            print_kv(&format!("ode {} y-hash", $name), hash_f64s(&y_flat));
        }};
    }
    run_explicit!("RK4", RK4);
    run_explicit!("RALS3", RALS3);
    run_explicit!("RALS4", RALS4);
    run_explicit!("RK5", RK5);
    // embedded methods that drive the step() arms touched by clippy.
    // Use a simple decay problem so every method converges within the
    // step-iteration budget.
    struct Decay;
    impl ODEProblem for Decay {
        fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) -> anyhow::Result<()> {
            dy[0] = -y[0];
            Ok(())
        }
    }
    macro_rules! run_embedded {
        ($name:expr, $integ:expr) => {{
            let (t, y) = BasicODESolver::new($integ)
                .solve(&Decay, (0.0, 2.0), 0.05, &[1.0])
                .unwrap();
            let y_flat: Vec<f64> = y.iter().flatten().copied().collect();
            print_kv(&format!("ode {} t-len", $name), t.len());
            print_kv(&format!("ode {} y-hash", $name), hash_f64s(&y_flat));
        }};
    }
    run_embedded!("RKF45", RKF45::new(1e-6, 0.9, 1e-6, 1e-1, 200));
    run_embedded!("DP45", DP45::new(1e-6, 0.9, 1e-6, 1e-1, 200));
    run_embedded!("TSIT45", TSIT45::new(1e-6, 0.9, 1e-6, 1e-1, 200));
    run_embedded!("RKF78", RKF78::new(1e-6, 0.9, 1e-6, 1e-1, 200));
    run_embedded!("BS23", BS23::new(1e-6, 0.9, 1e-6, 1e-1, 200));
    run_embedded!("GL4", GL4::new(ImplicitSolver::FixedPoint, 1e-6, 100));
}

fn check_distributions() {
    let mut rng = SmallRng::seed_from_u64(0xC0FFEE);
    // dist.rs:572 sample_with_rng arm
    let bern = Bernoulli(0.37);
    let s = bern.sample_with_rng(&mut rng, 100);
    print_kv("Bernoulli(0.37).sample[100]", hash_f64s(&s));
    // dist.rs:646 sample_with_rng arm (Uniform)
    let unif = Uniform(-2.5, 3.25);
    let s = unif.sample_with_rng(&mut rng, 100);
    print_kv("Uniform(-2.5,3.25).sample[100]", hash_f64s(&s));
}

fn check_matrix() {
    // py_matrix (api.rs:54-55)
    let m = py_matrix(vec![
        vec![1.0_f64, 2.0, 3.0, 4.0],
        vec![5.0, 6.0, 7.0, 8.0],
        vec![9.0, 10.0, 11.0, 12.0],
    ]);
    print_kv("py_matrix 3x4 hash", hash_matrix(&m));
    // change_shape (matrix.rs:897, 905)
    let m_col = m.change_shape();
    print_kv("change_shape Row->Col hash", hash_matrix(&m_col));
    let m_back = m_col.change_shape();
    print_kv("change_shape Col->Row hash", hash_matrix(&m_back));
    // diag (matrix.rs:1105)
    let sq = ml_matrix("1 2 3; 4 5 6; 7 8 9");
    let d = sq.diag();
    print_kv("diag 3x3", hash_f64s(&d));
    // to_vec (matrix.rs:1184)
    let v_of_v = sq.to_vec();
    let flat: Vec<f64> = v_of_v.into_iter().flatten().collect();
    print_kv("to_vec hash", hash_f64s(&flat));
    // col_reduce / row_reduce (matrix.rs:2863, 2874)
    let cr = sq.col_reduce(|c| c.iter().sum::<f64>());
    let rr = sq.row_reduce(|c| c.iter().sum::<f64>());
    print_kv("col_reduce sum", hash_f64s(&cr));
    print_kv("row_reduce sum", hash_f64s(&rr));
    // is_symmetric (matrix.rs:3660, 3683)
    let sym = ml_matrix("1 2 3; 2 5 6; 3 6 9");
    let nonsym = ml_matrix("1 2 3; 4 5 6; 7 8 9");
    print_kv("is_symmetric sym", sym.is_symmetric());
    print_kv("is_symmetric nonsym", nonsym.is_symmetric());
    // matrix.rs:764-765 - inspect; this site is inside `py_matrix`
    // creation path already exercised above. Cover additionally with a
    // non-square shape so both row/col loops fire.
    let rect = py_matrix(vec![vec![1.0_f64, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]]);
    print_kv("py_matrix 3x2 hash", hash_matrix(&rect));
}

fn check_polynomial() {
    let p = poly(vec![1.0, -2.0, 0.5, 3.0]); // 1 - 2x + 0.5x^2 + 3x^3
    print_kv("poly.derivative", hash_f64s(&p.derivative().coef));
    print_kv("poly.integral", hash_f64s(&p.integral().coef));
    for &x in &[-1.0_f64, -0.25, 0.0, 0.5, 1.5, 2.7] {
        print_kv(&format!("poly.eval({x})"), format!("{:.17e}", p.eval(x)));
    }
}

fn check_sparse() {
    // SPMatrix transpose (sparse.rs:97, 102)
    let m = ml_matrix("1 0 2 0; 0 3 0 4; 5 0 0 6; 0 7 8 0");
    let sp: SPMatrix = SPMatrix::from_dense(&m);
    let spt: SPMatrix = sp.transpose();
    let recovered = spt.to_dense();
    print_kv("sparse transpose hash", hash_matrix(&recovered));
    // Also exercise a non-square shape and verify double-transpose
    // returns the original.
    let m2 = ml_matrix("1 0 0; 0 2 3; 4 0 5; 0 6 0");
    let sp2: SPMatrix = SPMatrix::from_dense(&m2);
    let sp2tt = sp2.transpose().transpose().to_dense();
    print_kv("sparse round-trip hash", hash_matrix(&sp2tt));
}

fn check_useful() {
    // useful::max / useful::min (useful.rs:82, 100)
    let v = vec![3.2_f64, -1.0, 4.4, -7.0, 0.0, 10.5, -10.5, 2.2];
    print_kv("max(vec)", format!("{:.17e}", max(v.clone())));
    print_kv("min(vec)", format!("{:.17e}", min(v.clone())));
    let v2 = vec![1.0_f64, 2.0]; // l == 2 edge case
    print_kv("max([1,2])", format!("{:.17e}", max(v2.clone())));
    print_kv("min([1,2])", format!("{:.17e}", min(v2.clone())));
}

fn main() {
    println!("== integral ==");
    check_integral();
    println!();
    println!("== chebyshev/lagrange ==");
    check_chebyshev_lagrange();
    println!();
    println!("== lanczos ==");
    check_lanczos();
    println!();
    println!("== ode ==");
    check_ode();
    println!();
    println!("== distributions ==");
    check_distributions();
    println!();
    println!("== matrix ==");
    check_matrix();
    println!();
    println!("== polynomial ==");
    check_polynomial();
    println!();
    println!("== sparse ==");
    check_sparse();
    println!();
    println!("== useful ==");
    check_useful();
}
