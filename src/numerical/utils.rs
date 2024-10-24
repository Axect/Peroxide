use crate::structure::matrix::*;
use crate::traits::matrix::MatrixTrait;
use crate::structure::ad::*;
use crate::structure::ad::AD::*;
use crate::util::non_macro::{cat, zeros};

/// Jacobian Matrix
///
/// # Description
/// : Exact jacobian matrix using Automatic Differenitation
///
/// # Type
/// `(F, &Vec<f64>) -> Matrix where F: Fn(&Vec<AD>) -> Vec<AD>`
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
///     // r[0]    1   -1
///     // r[1]    1    2
/// }
///
/// fn f(xs: &Vec<AD>) -> Vec<AD> {
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
pub fn jacobian<F: Fn(&Vec<AD>) -> Vec<AD>>(f: F, x: &Vec<f64>) -> Matrix {
    let l = x.len();
    let mut x_ad: Vec<AD> = x.iter().map(|&x| AD1(x, 0f64)).collect();
    let l2 = f(&x_ad).len();

    let mut J = zeros(l2, l);

    for i in 0 .. l {
        x_ad[i][1] = 1f64;
        let slopes: Vec<f64> = f(&x_ad).iter().map(|ad| ad.dx()).collect();
        J.subs_col(i, &slopes);
        x_ad[i][1] = 0f64;
    }
    J
}

///// Hessian Matrix
//#[allow(non_snake_case)]
//pub fn hessian<F: Fn(&Vec<AD>) -> AD>(f: F, x: &Vec<f64>) -> Matrix {
//    let l = x.len();
//    let mut x_ad: Vec<AD> = x.iter().map(|&x| AD2(x, 0f64, 0f64)).collect();
//
//    let mut H = zeros(l, l);
//
//    for i in 0 .. l {
//        for j in 0 .. l {
//        }
//    }
//
//    unimplemented!()
//}

//#[allow(non_snake_case)]
//pub fn jacobian_ad<F: Fn(&Vec<AD>) -> Vec<AD>>(f: F, x: &Vec<AD>) -> Vec<Vec<AD>> {
//    let l = x.len();
//    let mut x_ad: Vec<AD> = x.clone().into_iter().map(|mut t| {
//        t.iter_mut().skip(1).for_each(|k| *k = 0f64);
//        t
//    }).collect();
//    let l2 = f(&x_ad).len();
//
//    let mut JT: Vec<Vec<AD>> = vec![vec![AD0(0f64); l2]; l];
//
//    for i in 0 .. l {
//        x_ad[i][1] = 1f64;
//        let ads = f(&x_ad);
//        JT[i] = ads;
//        x_ad[i][1] = 0f64;
//    }
//    JT
//}


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

    let a = cat(0f64, &a_input);
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
