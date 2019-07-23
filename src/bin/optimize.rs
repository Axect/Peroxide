extern crate peroxide;
use peroxide::*;

fn main() {
    // Make test data set
    let noise = Normal(0., 0.5);
    let p_true: Vec<Number> = NumberVector::from_f64_vec(vec![20f64, 10f64, 1f64, 50f64]);
    let real = f(p_true.clone()).to_f64_vec();
    let y = matrix(zip_with(|x, y| x + y, &real, &noise.sample(100)), 100, 1, Col);

    // Initial Value
    let p_init = vec![5f64, 2f64, 0.2f64, 10f64];
    let j_init = jacobian(f, p_init.clone());

    // Gradient Descent
    let mut p_gd = matrix(p_init.clone(), 4, 1, Col);
    let mut j_gd = j_init.clone();
    let mut y_hat_gd = matrix(f(NumberVector::from_f64_vec(p_init.clone())).to_f64_vec(), 100, 1, Col);

    for i in 0 .. 30 {
        let h = 0.001 / y.var()[0] * j_gd.t() * (&y - &y_hat_gd);
        p_gd = &p_gd + &h;
        j_gd = jacobian(f, p_gd.data.clone());
        y_hat_gd = matrix(f(NumberVector::from_f64_vec(p_gd.data.clone())).to_f64_vec(), 100, 1, Col);
    }

    p_gd.print();
}

// Non autonomous case
fn f(p: Vec<Number>) -> Vec<Number> {
    (0 .. 100).map(|t| Number::from_f64(t as f64)).map(|t| p[0] * (-t / p[1]).exp() + p[2] * t * (-t / p[3]).exp()).collect()
}

