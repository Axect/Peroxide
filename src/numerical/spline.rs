use operation::extra_ops::PowOps; // For pow
#[allow(unused_imports)]
use structure::matrix::*;
#[allow(unused_imports)]
use structure::polynomial::*;
#[allow(unused_imports)]
use structure::vector::*;
#[allow(unused_imports)]
use util::non_macro::*;

/// Cubic Spline (Natural)
///
/// # Description
///
/// Implement algorithm of Natural cubic splines, Arne Morten Kvarving.
///
/// # Type
/// (Vector, Vector) -> Vec<Polynomial>
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let x = c!(0.9, 1.3, 1.9, 2.1);
/// let y = c!(1.3, 1.5, 1.85, 2.1);
///
/// let s = cubic_spline(x, y);
///
/// for i in 0 .. s.len() {
///     s[i].print();
/// }
///
/// // -0.2347x^3 + 0.6338x^2 - 0.0329x + 0.9873
/// // 0.9096x^3 - 3.8292x^2 + 5.7691x - 1.5268
/// // -2.2594x^3 + 14.2342x^2 - 28.5513x + 20.2094
/// ```
pub fn cubic_spline(node_x: Vector, node_y: Vector) -> Vec<Polynomial> {
    //! Pre calculated variables
    //! node_x: n+1
    //! node_y: n+1
    //! h     : n
    //! b     : n
    //! v     : n
    //! u     : n
    //! z     : n+1
    let n = node_x.len() - 1;
    assert_eq!(n, node_y.len() - 1);

    // Pre-calculations
    let mut h = vec![0f64; n];
    let mut b = vec![0f64; n];
    let mut v = vec![0f64; n];
    let mut u = vec![0f64; n];
    for i in 0..n {
        if i == 0 {
            h[i] = node_x[i + 1] - node_x[i];
            b[i] = (node_y[i + 1] - node_y[i]) / h[i];
        } else {
            h[i] = node_x[i + 1] - node_x[i];
            b[i] = (node_y[i + 1] - node_y[i]) / h[i];
            v[i] = 2. * (h[i] + h[i - 1]);
            u[i] = 6. * (b[i] - b[i - 1]);
        }
    }

    // Tri-diagonal matrix
    let mut m = matrix(vec![0f64; (n - 1) * (n - 1)], n - 1, n - 1, Col);
    for i in 0..n - 2 {
        m[(i, i)] = v[i + 1];
        m[(i + 1, i)] = h[i + 1];
        m[(i, i + 1)] = h[i + 1];
    }
    m[(n - 2, n - 2)] = v[n - 1];

    // Calculate z
    let z_inner = m.inv().unwrap() * Vec::from(&u[1..]).to_matrix();
    let mut z = vec![0f64];
    z.extend(&z_inner.data);
    z.push(0f64);

    // Declare empty spline
    let mut s: Vec<Polynomial> = Vec::new();

    // Main process
    for i in 0..n {
        // Memoization
        let t_i = node_x[i];
        let t_i1 = node_x[i + 1];
        let z_i = z[i];
        let z_i1 = z[i + 1];
        let h_i = h[i];
        let y_i = node_y[i];
        let y_i1 = node_y[i + 1];
        let temp1 = poly(vec![1f64, -t_i]);
        let temp2 = poly(vec![1f64, -t_i1]);

        let term1 = temp1.powi(3) * (z_i1 / (6f64 * h_i));
        let term2 = temp2.powi(3) * (-z_i / (6f64 * h_i));
        let term3 = temp1 * (y_i1 / h_i - z_i1 * h_i / 6.);
        let term4 = temp2 * (-y_i / h_i + h_i * z_i / 6.0);

        s.push(term1 + term2 + term3 + term4);
    }
    return s;
}
