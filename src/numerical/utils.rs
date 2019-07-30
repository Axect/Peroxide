use structure::dual::*;
use structure::matrix::*;
use util::non_macro::{cat, zeros};
use ::{Number, NumberVector};

/// Jacobian Matrix
///
/// # Description
/// : Exact jacobian matrix using Automatic Differenitation
///
/// # Type
/// (Vector, F) -> Matrix where F: Fn(Vec<Dual>) -> Vec<Dual>
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let x = c!(1, 1);
/// let j = jacobian(f, x);
/// j.print();
///
/// //      c[0] c[1]
/// // r[0]    1    1
/// // r[1]   -1    2
///
/// fn f(xs: Vec<Number>) -> Vec<Number> {
///     let x = xs[0];
///     let y = xs[1];
///
///     vec![
///        x - y,
///        x + 2.*y,
///    ]
/// }
/// ```
#[allow(non_snake_case)]
pub fn jacobian<F>(g: F, x: Vec<f64>) -> Matrix
where
    F: Fn(Vec<Number>) -> Vec<Number>,
{
    let l = x.len();

    let f = |x: Vec<Dual>| g(NumberVector::from_dual_vec(x)).to_dual_vec();

    let x_var: Vec<Dual> = merge_dual(x.clone(), vec![1f64; l]);
    let x_const = x.clone().conv_dual();

    let l2 = f(x_const.clone()).len();

    let mut J = zeros(l2, l);

    let mut x_temp = x_const.clone();

    for i in 0..l {
        x_temp[i] = x_var[i];
        let dual_temp = f(x_temp.clone());
        let slope_temp = dual_temp.slopes();
        for j in 0..l2 {
            J[(j, i)] = slope_temp[j];
        }
        x_temp = x_const.clone();
    }
    J
}

/// TriDiagonal Matrix Algorithm (TDMA)
///
/// # Description
///
/// Solve tri-diagonal matrix system efficiently (O(n))
/// ```bash
/// |b0 c0         | |x0|   |y0|
/// |a1 b1 c1      | |x1|   |y1|
/// |   a2 b2 c2   | |x2| = |y2|
/// |      ...     | |..|   |..|
/// |         am bm| |xm|   |ym|
/// ```
///
/// # Caution
///
/// You should apply boundary condition yourself
pub fn tdma(a_input: Vec<f64>, b_input: Vec<f64>, c_input: Vec<f64>, y_input: Vec<f64>) -> Matrix {
    let n = b_input.len();
    assert_eq!(a_input.len(), n - 1);
    assert_eq!(c_input.len(), n - 1);
    assert_eq!(y_input.len(), n);

    let a = cat(0f64, a_input.clone());
    let mut b = b_input.clone();
    let c = {
        let mut c_temp = c_input.clone();
        c_temp.push(0f64);
        c_temp.clone()
    };
    let mut y = y_input.clone();

    // Forward substitution
    let mut w = vec![0f64; n];
    for i in 1..n {
        w[i] = a[i] / b[i - 1];
        b[i] = b[i] - w[i] * c[i - 1];
        y[i] = y[i] - w[i] * y[i - 1];
    }

    // Backward substitution
    let mut x = vec![0f64; n];
    x[n - 1] = y[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (y[i] - c[i] * x[i + 1]) / b[i];
    }
    x.to_matrix()
}
