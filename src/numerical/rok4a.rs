#[allow(unused_imports)]
use structure::matrix::*;
use structure::dual::*;
use structure::vector::*;
use numerical::utils::jacobian;
use util::non_macro::{zeros, concat, eye, cbind};
use util::print::*;

#[allow(non_snake_case)]
pub fn one_step_rok4a<F>(xs: Vec<f64>, f:F, step: f64) -> Vec<f64>
    where F: Fn(Vec<Dual>) -> Vec<Dual> + Copy
{
    let ALPHA_MAT = ALPHA_MAT();
    let GAMMA_MAT = GAMMA_MAT();

    let n = xs.len() - 1;
    let m = 4;
    let s = 4;
    let t = xs[0];
    let h = step;

    let ys = xs.clone()
        .into_iter()
        .skip(1)
        .collect::<Vec<f64>>();

    let (H, V, w) = modified_arnoldi(xs.clone(), f);
    // F1 = n x 1
    let F1 = matrix(
        f(xs.conv_dual()).values(),
        n,
        1,
        Col
    );
    // phi = M x 1
    let phi1 = V.t() % F1.clone() + w.clone();
    let mult1 = (eye(m) - h*GAMMA*H.clone()).inv()
        .unwrap();
    // lambda = M x 1
    let lambda1 = mult1 % (h*phi1.clone());
    // k0 = n x 1
    let k0 = V.clone() % lambda1.clone() + h * (F1 - V.clone() % phi1);

    let mut k = k0;
    let mut lambda = lambda1;

    for i in 1 .. s {
        let new_t = t + ALPHA_VEC[i]*step;
        let mut s_vec = vec![0f64; n];
        for j in 0 .. i {
            let aij = ALPHA_MAT[(i, j)];
            let target = k.col(j).fmap(|x| x * aij);
            s_vec = s_vec.add(&target);
        }
        let new_y = ys.add(&s_vec);

        // F: n x 1
        let F = matrix(
            f(concat(vec![new_t], new_y).conv_dual()).values(),
            n,
            1,
            Col
        );
        let phi = V.t() % F.clone() + w.clone();
        let mult1 = (eye(m) - h*GAMMA*H.clone()).inv()
            .unwrap();
        let mut s_vec2 = vec![0f64; m];
        for j in 0 .. i {
            let gij = GAMMA_MAT[(i, j)];
            let target = lambda.col(j).fmap(|x| x * gij);
            s_vec2 = s_vec2.add(&target);
        }
        let mult2 = h * (phi.clone() + (H.clone() % s_vec2));
        let lambda_i = mult1 % mult2;
        let k_i = V.clone() % lambda_i.clone() + h * (F - V.clone() % phi);

        k = cbind(k, k_i);
        lambda = cbind(lambda, lambda_i);
    }

    let new_t = t + h;
    let mut s_vec = vec![0f64; n];
    for i in 0 .. s {
        let target = k.col(i).fmap(|x| x * B.clone()[i]);
        s_vec = s_vec.add(&target);
    }

    let new_ys = ys.add(&s_vec);
    concat(vec![new_t], new_ys)
}


#[allow(non_snake_case)]
pub fn modified_arnoldi<F>(xs: Vec<f64>, f:F) -> (Matrix, Matrix, Matrix)
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
    let mut w: f64 = 1f64 / beta;
    let mut v: Vec<f64> = f_n.fmap(|x| x / beta);

    let m: usize = 4;

    let mut ws: Vec<f64> = vec![0f64; m];
    ws[0] = w;

    // J = ∂f/∂(t,x) = n x (n+1)
    // J_n = ∂f/∂y = n x n
    let J_temp = jacobian(xs.clone(), f);
    let mut J = zeros(n,n);

    for i in 0 .. n {
        for j in 0 .. n {
            J[(i, j)] = J_temp[(i, j+1)];
        }
    }

    // Only differentiation with t
    let t_dual = dual(t, 1.);
    let ys_dual = ys.conv_dual();

    let xs_target = concat(vec![t_dual], ys_dual);

    let mut H = zeros(m, m);
    let mut V = zeros(n, m);

    let kappa: f64 = 0.25;

    for i in 0 .. m {
        // Fill V
        for j in 0 .. n {
            V[(j, i)] = v[j];
        }
        let f_t: Vec<f64> = f(xs_target.clone()).slopes();
        // ζ = J_n v_i + f_t(t_n, y_n)ω_i
        let mut zeta: Vec<f64> = (J.clone() % v.clone()).col(0).add(&(f_t.fmap(|x| x * w)));
        // ξ = 0, τ = |ζ|
        let mut xi = 0f64;
        let tau = zeta.norm();

        for j in 0 .. (i+1) {
            let vj = V.col(j);
            // H_{ji} = <ζ,v_j> + ξω_j
            H[(j, i)] = zeta.dot(&vj) + xi * ws[j];
            let hji = H[(j, i)];

            // ζ = ζ - H_{j,i}v_j, ξ = ξ - H_{j,i}ω_j
            let dzeta = vj.fmap(|v| hji * v);
            let dxi = hji * ws[j];
            zeta = zeta.sub(&dzeta);
            xi = xi - dxi;
        }

        let det_vec = concat(zeta.clone(), vec![xi]);
        let det_norm = det_vec.norm() / tau;
        
        // |[ζ^T ξ]^T| / τ <= κ
        if det_norm <= kappa {
            for j in 0 .. (i+1) {
                let vj = V.col(j);
                // ρ = <ζ, v_j> + ξω_j
                let rho = zeta.dot(&vj) + xi * ws[j];

                // ζ = ζ - ρv_j, ξ = ξ - ρω_j
                let dzeta2 = vj.fmap(|v| rho * v);
                let dxi2 = rho * ws[j];
                zeta = zeta.sub(&dzeta2);
                xi = xi - dxi2;
                H[(j, i)] = H[(j, i)] + rho;
            }
        }

        let det_vec2 = concat(zeta.clone(), vec![xi]);
        let det_norm2 = det_vec2.norm();
        let hii = det_norm2;

        v = zeta.fmap(|z| z / hii);
        w = xi / hii;

        if i < (m-1) {
            H[(i+1, i)] = det_norm2;
            ws[i+1] = w;
        }
    }

    (H, V, matrix(ws, m, 1, Col))
}

const GAMMA: f64 = 0.572816062482135;

const ALPHA21: f64 = 1.;
const ALPHA31: f64 = 0.10845300169319391758;
const ALPHA32: f64 = 0.39154699830680608241;
const ALPHA41: f64 = 0.43453047756004477624;
const ALPHA42: f64 = 0.14484349252001492541;
const ALPHA43: f64 = -0.07937397008005970166;

const ALPHA2: f64 = 1.;
const ALPHA3: f64 = 0.5;
const ALPHA4: f64 = 0.5;

const ALPHA_VEC: [f64;4] = [0., 1., 0.5, 0.5];

#[allow(non_snake_case)]
pub fn ALPHA_MAT() -> Matrix {
    Matrix {
        data: vec![
            0.,      0.,      0.,      0.,
            ALPHA21, 0.,      0.,      0.,
            ALPHA31, ALPHA32, 0.,      0.,
            ALPHA41, ALPHA42, ALPHA43, 0.,
        ],
        row: 4,
        col: 4,
        shape: Row,
    }
}

const GAMMA21: f64 = -1.91153192976055097824;
const GAMMA31: f64 = 0.32881824061153522156 ;
const GAMMA32: f64 = 0.;
const GAMMA41: f64 = 0.03303644239795811290;
const GAMMA42: f64 = -0.24375152376108235312;
const GAMMA43: f64 = -0.17062602991994029834;

#[allow(non_snake_case)]
pub fn GAMMA_MAT() -> Matrix {
    Matrix {
        data: vec![
            0., 0., 0., 0.,
            GAMMA21, 0., 0., 0.,
            GAMMA31, GAMMA32, 0., 0.,
            GAMMA41, GAMMA42, GAMMA43, 0.,
        ],
        row: 4,
        col: 4,
        shape: Row,
    }
}

const B1: f64 = 1f64/6f64;
const B2: f64 = 1f64/6f64;
const B3: f64 = 0.;
const B4: f64 = 2f64/3f64;

const B: [f64;4] = [B1, B2, B3, B4];
