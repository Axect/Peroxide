extern crate peroxide;
use peroxide::*;

const SIZE: usize = 100;

fn main() {
    let size = SIZE;

    // Make test data set
    let noise = Normal(0., 0.5);
    let p_true: Vec<Number> = NumberVector::from_f64_vec(vec![20f64, 10f64, 1f64, 50f64]);
    let real = f(p_true.clone()).to_f64_vec();
    let y = matrix(zip_with(|x, y| x + y, &real, &noise.sample(size)), size, 1, Col);
    y.print();
    real.print();

    // Initial Value
    let p_init = vec![5f64, 2f64, 0.2f64, 10f64];
    let j_init = jacobian(f, p_init.clone());

    // Gradient Descent
    let mut p_gd = matrix(p_init.clone(), 4, 1, Col);
    let mut j_gd = j_init.clone();
    let mut y_hat_gd = matrix(f(NumberVector::from_f64_vec(p_init.clone())).to_f64_vec(), size, 1, Col);

    for i in 0 .. 30 {
        let h = 0.02 * j_gd.t() * (&y - &y_hat_gd);
        p_gd = &p_gd + &h;
        j_gd = jacobian(f, p_gd.data.clone());
        y_hat_gd = matrix(f(NumberVector::from_f64_vec(p_gd.data.clone())).to_f64_vec(), size, 1, Col);
    }

    p_gd.print();

    // Gauss_Newton
    let mut p_gn = matrix(p_init.clone(), 4, 1, Col);
    let mut y_hat_gn = matrix(f(NumberVector::from_f64_vec(p_init.clone())).to_f64_vec(), size, 1, Col);
    let mut j_gn = j_init.clone();
    for i in 0 .. 10 {
        let h_gn: Matrix;
        match j_gn.pseudo_inv() {
            Some(W) => h_gn = W * (&y - &y_hat_gn),
            None => break,
        }
        p_gn = &p_gn + &h_gn;
        j_gn = jacobian(f, p_gn.data.clone());
        y_hat_gn = matrix(f(NumberVector::from_f64_vec(p_gn.data.clone())).to_f64_vec(), size, 1, Col);
    }

    p_gn.print();

    // Levenberg-Marquardt
    let (lambda_0, eps1, eps2) = (1e-2, 1e-6, 1e-6);
    let mut p_lm = p_init.to_matrix();
    let mut y_hat_lm = f(NumberVector::from_f64_vec(p_lm.data.clone())).to_f64_vec().to_matrix();
    let mut j_lm = j_init.clone();
    let mut A = &j_lm.t() * &j_lm;
    let mut lambda = lambda_0 * max(A.diag());
    let mut chi2 = ((&y - &y_hat_lm).t() * (&y - &y_hat_lm))[(0,0)];
    let mut nu = 2f64;
    
    for i in 0 .. 30 {
        let h_lm: Matrix;
        let mut A_diag = eye(A.row);
        let A_ref = A.diag();
        for i in 0 .. A.row {
            A_diag[(i,i)] = A_ref[i];
        }
        match (A.clone() + lambda * A_diag.clone()).inv() {
            Some(B) => h_lm = B * j_lm.t() * (&y - &y_hat_lm),
            None => break,
        }

        let p_lm_temp = &p_lm + &h_lm;
        let j_lm_temp = jacobian(f, p_lm.data.clone());
        let A_temp = &j_lm.t() * &j_lm;
        let y_hat_temp = f(NumberVector::from_f64_vec(p_lm_temp.data.clone())).to_f64_vec().to_matrix();
        let chi2_temp = ((&y - &y_hat_temp).t() * (&y - &y_hat_temp))[(0,0)];

        let rho = (chi2 - chi2_temp) / (h_lm.t() * (lambda * A_diag * h_lm.clone() + j_lm.t() * (&y - &y_hat_lm)))[(0,0)];

        rho.print();

        if rho > 0f64 {
            p_lm = p_lm_temp;
            j_lm = j_lm_temp;
            A = A_temp;
            y_hat_lm = y_hat_temp;
            chi2 = chi2_temp;

            lambda = lambda * max(c!(1f64/3f64, 1f64 - (2f64*rho - 1f64).powi(3)));
            nu = 2f64;
        } else {
            lambda = lambda * nu;
            nu *= 2f64;
        }
    }

    p_lm.print();

    // Plot
    let p_x = seq(0, (SIZE-1) as i32, 1);
    let p_y = y.data;
    let f_y = f(NumberVector::from_f64_vec(p_lm.data)).to_f64_vec();

    let mut plt = Plot2D::new();
    plt.set_domain(p_x)
        .insert_image(p_y)
        .insert_image(f_y)
        .set_legends(vec!["real", "fit"])
        .set_fig_size((10, 6))
        .set_dpi(300)
        .set_path("example_data/lm.png")
        .set_title("Levenberg-Marquardt Algorithm")
        .set_xlabel("$x$")
        .set_ylabel("$y$")
        .savefig();
}

// Non autonomous case
fn f(p: Vec<Number>) -> Vec<Number> {
    (0 .. SIZE).map(|t| Number::from_f64(t as f64)).map(|t| p[0] * (-t / p[1]).exp() + p[2] * t * (-t / p[3]).exp()).collect()
}

