use structure::matrix::*;
use structure::vector::*;
use structure::dual::*;
use numerical::utils::jacobian;
use util::non_macro::concat;

/// Value of 3f64.sqrt()
const SQRT3: f64 = 1.7320508075688772;

/// Butcher tableau for Gauss_legendre 4th order
const GL4: [[f64; 3]; 2] = [
    [0.5 - SQRT3/6f64, 0.25, 0.25 - SQRT3/6f64],
    [0.5 + SQRT3/6f64, 0.25 + SQRT3/6f64, 0.25]
];

pub fn newton_iter<F>(f: F, t: f64, h: f64, y: Vec<f64>, kl: Vec<f64>) -> Vec<Dual>
    where F: Fn(Dual, Vec<Dual>) -> Vec<Dual> + Copy,
{
    let t1 = dual(t + GL4[0][0] * h, 0.);
    let t2 = dual(t + GL4[1][0] * h, 0.);
    let yn = y.conv_dual();
    let n = y.len();

    let g = |k: Vec<Dual>| -> Vec<Dual> {
        let k1 = k.take(n);
        let k2 = k.skip(n);
        concat(
            f(t1, yn.add(&k1.fmap(|x| x * GL4[0][1] * h).add(&k2.fmap(|x| x * GL4[0][2] * h)))),
            f(t2, yn.add(&k1.fmap(|x| x * GL4[1][1] * h).add(&k2.fmap(|x| x * GL4[1][2] * h))))
        )
    };

    unimplemented!()
}