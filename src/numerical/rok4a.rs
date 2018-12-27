#[allow(unused_imports)]
use structure::matrix::*;
use structure::dual::*;
use structure::vector::*;
use numerical::utils::jacobian;
use util::non_macro::{zeros, concat};

#[allow(non_snake_case)]
pub fn modified_arnoldi<F>(xs: Vec<f64>, f:F)
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    // Non-autonomous Equation
    // y' = f(xs) = f(t, ys)
    
    let n = xs.len() - 1;

    // (t, ys) = xs
    let t: f64 = xs[0];
    let ys: Vec<f64> = xs.clone()
        .into_iter()
        .skip(1)
        .collect::<Vec<f64>>();

    // f(xs) - Vec<f64>
    let f_n: Vec<f64> = f(xs.conv_dual()).conv_dual();

    // β = |[f_n^T 1]^T|
    // ω = 1 / β
    // v = f_n / β
    let beta: f64 = concat(f_n.clone(), vec![1f64]).norm();
    let mut omega: f64 = 1f64 / beta;
    let mut v: Vec<f64> = f_n.fmap(|x| x / beta);

    let m: usize = 4;

    // J = ∂f/∂(t,x) = n x (n+1)
    // J_n = ∂f/∂y = n x n
    let J_temp = jacobian(ys.clone(), f);
    let mut J = zeros(n,n);

    for i in 0 .. n {
        for j in 0 .. n {
            J[(i, j)] = J_temp[(i, j+1)];
        }
    }
    
    // Only differentiation with t
    let t_dual = dual(t, 1.);
    let ys_dual = ys.conv_dual();

    let mut xs_target = concat(vec![t_dual], ys_dual);

    let mut H = zeros(m, m);
    let mut V = zeros(m ,m);

    for i in 0 .. m {
        // Fill V
        for j in 0 .. m {
            V[(j, i)] = v[j];
        }
        let f_t: Vec<f64> = f(xs_target.clone()).slopes();
        let zeta: Vec<f64> = (J.clone() % v.clone()).col(0).add(&(f_t.fmap(|x| x * omega)));
        let xi = 0f64;
        let tau = zeta.norm();

        for j in 0 .. (i+1) {
            H[(j, i)] = zeta.dot(&V.col(j));
        }
    }

    unimplemented!()
}






const γ: f64 = 0.572816062482135;

const α21: f64 = 1.;
const α31: f64 = 0.10845300169319391758;
const α32: f64 = 0.39154699830680608241;
const α41: f64 = 0.43453047756004477624;
const α42: f64 = 0.14484349252001492541;
const α43: f64 = -0.07937397008005970166;

const γ21: f64 = -1.91153192976055097824;
const γ31: f64 = 0.32881824061153522156 ;
const γ32: f64 = 0.;
const γ41: f64 = 0.03303644239795811290;
const γ42: f64 = -0.24375152376108235312;
const γ43: f64 = -0.17062602991994029834;

const b1: f64 = 1f64/6f64;
const b2: f64 = 1f64/6f64;
const b3: f64 = 0.;
const b4: f64 = 2f64/3f64;
