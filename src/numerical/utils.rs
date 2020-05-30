use crate::structure::dual::*;
use crate::structure::matrix::*;
use crate::util::non_macro::{cat, zeros};
use crate::traits::num::{Real, Number, NumberVector};

/// Jacobian Matrix
///
/// # Description
/// : Exact jacobian matrix using Automatic Differenitation
///
/// # Type
/// (Vector, F) -> Matrix where F: Fn(Vec<Number>) -> Vec<Number>
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// 
/// fn main() {
///     let x = c!(1, 1);
///     let j = jacobian(f, &x);
///     j.print();
///
///     //      c[0] c[1]
///     // r[0]    1    1
///     // r[1]   -1    2
/// }
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
pub fn jacobian<F>(g: F, x: &Vec<f64>) -> Matrix
where
    F: Fn(Vec<Number>) -> Vec<Number>,
{
    let f = |x: &Vec<Dual>| g(NumberVector::from_dual_vec(x.clone())).to_dual_vec();
    jacobian_real(Box::new(f), x)
}

/// Jacobian Matrix for Real input
///
/// # Description
/// : Exact jacobian matrix using Automatic Differenitation
///
/// # Type
/// (Box<F>, &Vec<T>) -> Matrix where F: Fn(&Vec<Dual>) -> Vec<Dual>, T: Real
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// 
/// fn main() {
///     let x = c!(1, 1);
///     let j = jacobian_real(Box::new(f), &x);
///     j.print();
///
///     //      c[0] c[1]
///     // r[0]    1    1
///     // r[1]   -1    2
/// }
///
/// fn f(xs: &Vec<Dual>) -> Vec<Dual> {
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
pub fn jacobian_real<F, T>(f: Box<F>, x: &Vec<T>) -> Matrix
where
    T: Real,
    F: Fn(&Vec<Dual>) -> Vec<Dual>,
{
    let l = x.len();
    let mut x_dual: Vec<Dual> = x
        .clone()
        .into_iter()
        .map(|t| dual(t.to_f64(), 0f64))
        .collect();
    let l2 = f(&x_dual).len();

    let mut J = zeros(l2, l);

    for i in 0..l {
        x_dual[i].set_slope(1f64);
        let slopes = f(&x_dual).slopes();
        J.subs_col(i, &slopes);
        x_dual[i].set_slope(0f64);
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
    x.into()
}
