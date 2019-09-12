extern crate peroxide;
use peroxide::*;

#[test]
fn test_jacobian() {
    let x = c!(1, 0);
    let j = jacobian(f, &x);
    assert_eq!(j, ml_matrix("0 1; 5 1"));
}

/// Test function
///
/// # Function
/// $$\begin{pmatrix} x^2 y \\\ 5x+ \sin y \end{pmatrix}$$
///
/// # Jacobian
/// $$\begin{pmatrix} 2xy & x^2 \\\ 5 & \cos y \end{pmatrix}$$
fn f(xs: Vec<Number>) -> Vec<Number> {
    let x = xs[0];
    let y = xs[1];
    vec![x.powi(2) * y, 5f64 * x + y.sin()]
}
