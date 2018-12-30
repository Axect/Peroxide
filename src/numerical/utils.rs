use structure::matrix::*;
use structure::dual::*;
use util::non_macro::zeros;

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
/// let j = jacobian(x, f);
/// j.print();
/// 
/// //      c[0] c[1]
/// // r[0]    1    1
/// // r[1]   -1    2
/// 
/// fn f(xs: Vec<Dual>) -> Vec<Dual> {
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
pub fn jacobian<F>(x: Vec<f64>, f: F) -> Matrix
    where F: Fn(Vec<Dual>) -> Vec<Dual>
{
    let l = x.len();

    let x_var: Vec<Dual> = merge_dual(x.clone(), vec![1f64; l]);
    let x_const = x.clone().conv_dual();

    let l2 = f(x_const.clone()).len();

    let mut J = zeros(l2, l);

    let mut x_temp = x_const.clone();

    for i in 0 .. l {
        x_temp[i] = x_var[i];
        let dual_temp = f(x_temp.clone());
        let slope_temp = dual_temp.slopes();
        for j in 0 .. l2 {
            J[(j, i)] = slope_temp[j];
        }
        x_temp = x_const.clone();
    }
    J
}

/// Non Autonomous Jacobian
///
/// # Description
/// : For f(t, y), `J_{ij}=∂f_i/∂y_j`
#[allow(non_snake_case)]
pub fn non_auto_jacobian<F>(xs: Vec<f64>, f: F) -> Matrix
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    // n = len(ys) where xs = (t, ys)
    let n = xs.len() - 1;

    // original J = n x (n+1)
    let J_temp = jacobian(xs.clone(), f);

    // J for non autonomous
    // J_y = n x n
    let mut J = zeros(n, n);

    for i in 0 .. n {
        for j in 0 .. n {
            J[(i, j)] = J_temp[(i, j+1)];
        }
    }
    J
}
