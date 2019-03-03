extern crate peroxide;
use peroxide::*;

// x : n x L
// xb: n x (L+1)
// v : (L+1) x M
// a : n x M
// ab: n x (M+1)
// w : (M+1) x n
// wb: M x N
// y : n x N
// t : n x N
// dh: n x M
// do: n x N

fn main() {
    let v = weights_init(3, 2);
    let w = weights_init(3, 1);

    let x = ml_matrix("0 0; 0 1; 1 0; 1 1");
    let t = ml_matrix("0;1;1;0");

    let y = train(v, w, x, t, 0.25, 5000);
    y.print();
}

fn weights_init(m: usize, n: usize) -> Matrix {
    rand(m, n) * 2f64 - 1f64
}

fn sigmoid(x: f64) -> f64 {
    1f64 / (1f64 + (-x).exp())
}

fn forward(weights: Matrix, input_bias: Matrix) -> Matrix {
    let s = input_bias % weights;
    s.fmap(|x| sigmoid(x))
}

fn add_bias(input: Matrix, bias: f64) -> Matrix {
    let b = matrix(vec![bias; input.row], input.row, 1, Col);
    cbind(b, input)
}

fn hide_bias(weight: Matrix) -> Matrix {
    weight.skip(1, Row)
}

fn train(weights1: Matrix, weights2: Matrix, input: Matrix, answer: Matrix, eta: f64, times: usize) -> Matrix {
    let x = input;
    let mut v = weights1;
    let mut w = weights2;
    let t = answer;
    let xb = add_bias(x.clone(), -1f64);

    for i in 0 .. times {
        let a = forward(v.clone(), xb.clone());
        let ab = add_bias(a.clone(), -1f64);
        let y = forward(w.clone(), ab.clone());
        let err = (y.clone() - t.clone()).t() % (y.clone() - t.clone());
        let wb = hide_bias(w.clone());
        let delta_o = (y.clone() - t.clone()) * y.clone() * (1f64 - y.clone());
        let delta_h = (delta_o.clone() % wb.t()) * a.clone() * (1f64 - a.clone());

        w = w.clone() - eta * (ab.t() % delta_o);
        v = v.clone() - eta * (xb.t() % delta_h);
    }

    let a = forward(v, xb);
    let ab = add_bias(a, -1f64);
    let y = forward(w, ab);

    y
}