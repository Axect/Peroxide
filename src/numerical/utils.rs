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

    let mut J = zeros(l, l);

    let mut x_temp = x_const.clone();

    for i in 0 .. l {
        x_temp[i] = x_var[i];
        let dual_temp = f(x_temp.clone());
        let slope_temp = dual_temp.slopes();
        for j in 0 .. l {
            J[(j, i)] = slope_temp[j];
        }
        x_temp = x_const.clone();
    }
    J
}