extern crate peroxide;
use peroxide::*;

const N: usize = 50;

#[allow(non_snake_case)]
fn main() {
    let cos_45 = 2f64.sqrt() / 2f64;
    let R = matrix(vec![cos_45, -cos_45, cos_45, cos_45], 2, 2, Row);

    // Group 1
    let m1 = Coord { x: -1f64, y: 2f64 };
    let x1_x_temp = Normal(0f64, 1f64).sample(N);
    let x1_y_temp = Normal(0f64, 0.3f64).sample(N);
    let x1_temp = cbind(x1_x_temp.to_matrix(), x1_y_temp.to_matrix()).t();
    let x1_rot = R.clone() * x1_temp;
    let x1_x = x1_rot.row(0).fmap(|t| t + m1.x);
    let x1_y = x1_rot.row(1).fmap(|t| t + m1.y);
    let x1 = cbind(x1_x.to_matrix(), x1_y.to_matrix());
    x1.write_with_header("example_data/sample1.csv", vec!["x", "y"])
        .expect("Can't write file");

    // Group 2
    let m2 = Coord { x: 2f64, y: -2f64 };
    let x2_x_temp = Normal(0f64, 1f64).sample(N);
    let x2_y_temp = Normal(0f64, 0.3f64).sample(N);
    let x2_temp = cbind(x2_x_temp.to_matrix(), x2_y_temp.to_matrix()).t();
    let x2_rot = R.clone() * x2_temp;
    let x2_x = x2_rot.row(0).fmap(|t| t + m2.x);
    let x2_y = x2_rot.row(1).fmap(|t| t + m2.y);
    let x2 = cbind(x2_x.to_matrix(), x2_y.to_matrix());
    x2.write_with_header("example_data/sample2.csv", vec!["x", "y"])
        .expect("Can't write file");

    let X = rbind(x1.clone(), x2.clone());
    let X_bias = matrix(vec![1f64; 2 * N], 2 * N, 1, Col);
    let X_tilde = cbind(X_bias, X.clone());

    let mut T = zeros(2 * N, 2);
    T.subs_col(0, concat(vec![0f64; N], vec![1f64; N]));
    T.subs_col(1, concat(vec![1f64; N], vec![0f64; N]));

    x1.print();
    x2.print();
    X.print();

    X_tilde.print();

    let weight = X_tilde.pseudo_inv().unwrap() * T;

    weight.print();

    let line = line_function(weight.clone(), seq(-4, 4, 0.1));
    line.write_with_header("example_data/line.csv", vec!["x", "y"])
        .expect("Can't write file");

    let classifier = fisher_lda(weight);
    classifier(ml_matrix("-2; 2")).print();
    classifier(ml_matrix("2; -4")).print();
}

struct Coord {
    x: f64,
    y: f64,
}

#[allow(non_snake_case)]
fn fisher_lda(weight: Matrix) -> impl Fn(Matrix) -> Matrix {
    move |x: Matrix| {
        let x_tilde = add_bias(x, Row);
        weight.t() * x_tilde
    }
}

fn add_bias(x: Matrix, shape: Shape) -> Matrix {
    match shape {
        Col => {
            let n = x.row;
            let x_bias = matrix(vec![1f64; n], n, 1, Col);
            cbind(x_bias, x)
        }
        Row => {
            let n = x.col;
            let x_bias = matrix(vec![1f64; n], 1, n, Row);
            rbind(x_bias, x)
        }
    }
}

fn line_function(weight: Matrix, input: Vec<f64>) -> Matrix {
    let a_vec = weight.col(0);
    let ys = input.fmap(|x| -a_vec[1] / a_vec[2] * x + (0.5 - a_vec[0]) / a_vec[2]);
    cbind(input.to_matrix(), ys.to_matrix())
}
