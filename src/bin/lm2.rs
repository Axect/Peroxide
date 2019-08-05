extern crate peroxide;
use peroxide::*;

// =============================================================================
// Type Aliases
// =============================================================================
type Domain = Vec<f64>;
type Parameter = Vec<Number>;

fn main() {
    // =========================================================================
    // Domain Generation
    // =========================================================================
    let t = seq(0, 99, 1);
    let p_true: Vec<Number> = NumberVector::from_f64_vec(vec![20f64, 10f64, 1f64, 50f64]);
    let noise = Normal(0., 0.5);
    let y = &f(&t, p_true).to_f64_vec().to_matrix() + &noise.sample(100).to_matrix();

    // =========================================================================
    // Initialize
    // =========================================================================
    let mut p: Vec<Number> = NumberVector::from_f64_vec(vec![5f64, 2f64, 0.2, 10f64]);
    let mut p_mat = p.to_f64_vec().to_matrix();
    let mut y_hat = f(&t, p.clone()).to_f64_vec().to_matrix();
    let mut j = jacobian(|p| f(&t, p), p.to_f64_vec());
    let mut jtj = &j.t() * &j;
    let mut lambda = 1e-2 * max(jtj.diag());
    let mut nu = 2f64;

    // =========================================================================
    // Update
    // =========================================================================
    for _i in 0 .. 30 {
        let b = &jtj + &(lambda * jtj.to_diag());
        let h: Matrix;
        match b.inv() {
            None => break,
            Some(m) => {
                h = m * j.t() * (&y - &y_hat);
            }
        }
        let p_temp = &p_mat + &h;
        let p_temp_num: Vec<Number> = NumberVector::from_f64_vec(p_temp.col(0));
        let y_hat_temp = f(&t, p_temp_num.clone()).to_f64_vec().to_matrix();
        let rho = ((&y - &y_hat_temp).t() * (&y - &y_hat_temp))[(0,0)] / (&h.t() * &(lambda * (&jtj.to_diag() * &h) + j.t() * (&y - &y_hat)))[(0,0)];
        if rho > 0f64 {
            p_mat = p_temp;
            p = NumberVector::from_f64_vec(p_mat.col(0));
            y_hat = y_hat_temp;
            j = jacobian(|p| f(&t, p), p.to_f64_vec());
            jtj = &j.t() * &j;
            lambda = lambda * max(vec![1f64/3f64, 1f64 - (2f64*rho - 1f64).powi(3)]);
            nu = 2f64;
        } else {
            lambda = lambda * nu;
            nu *= 2f64;
        }
    }
    p_mat.print();
}

fn f(t: &Domain, p: Parameter) -> Parameter {
    t.clone().into_iter()
        .map(|x| p[0] * (-x / p[1]).exp() + p[2] * x * (-x / p[3]).exp())
        .collect()
}
