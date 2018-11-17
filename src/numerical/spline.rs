#[allow(unused_imports)]
use structure::matrix::*;
#[allow(unused_imports)]
use structure::vector::*;
#[allow(unused_imports)]
use structure::polynomial::*;
#[allow(unused_imports)]
use util::non_macro::*;

// pub fn monotone_cubic_spline(node_x: Vector, node_y: Vector) -> Polynomial {
//     let l = node_x.len();
//     assert_eq!(l, node_y.len());

//     if l == 0 {
//         return poly(c(&[0]));
//     } else if l == 1 {
//         return poly(node_y);
//     }



//     unimplemented!()
// }

fn secant(x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    (y2 - y1) / (x2 - x1)
}