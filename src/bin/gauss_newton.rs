extern crate peroxide;
use peroxide::*;

const S: [f64; 7] = [
    0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740
];

#[allow(non_snake_case)]
fn main() {
    let y = ml_matrix("0.05; 0.127; 0.094; 0.2122; 0.2729; 0.2665; 0.3317");

    let beta_init = vec!(0.9, 0.2);

    let mut beta = beta_init.to_matrix();
    let mut j = jacobian(f, beta_init);
    let mut y_hat = f(NumberVector::from_f64_vec(beta.data.clone()))
        .to_f64_vec()
        .to_matrix();

    for _i in 0 .. 10 {
        let h: Matrix;
        match j.pseudo_inv() {
            Some(W) => h = W * (&y - &y_hat),
            None => break,
        }
        beta = &beta + &h;
        j = jacobian(f, beta.data.clone());
        y_hat = f(NumberVector::from_f64_vec(beta.data.clone())).to_f64_vec().to_matrix();
    }

    beta.print();
}

fn f(beta: Vec<Number>) -> Vec<Number> {
    let s: Vec<Number> = NumberVector::from_f64_vec(S.to_vec());
    map(|x| beta[0] * x / (beta[1] + x), &s)
}