extern crate peroxide;
use peroxide::*;

fn main() {
    let mut j = jacobian(f, vec![0f64, 0f64, 0f64]);
    j.print();

    j = jacobian(f, vec![0f64, 1f64, 1f64]);
    j.print();

    j = jacobian(g, vec![1f64, 1f64, 0f64]);
    j.print();
    j.det().print();
}

fn f(v: Vec<Number>) -> Vec<Number> {
    vec![v[0], 5f64*v[2], 4f64*v[1].powi(2) - 2f64 * v[2], v[2] * v[0].sin()]
}

fn g(v: Vec<Number>) -> Vec<Number> {
    vec![5f64 * v[1], 4f64 * v[0].powi(2) - 2f64 * (v[1] * v[2]).sin(), v[1]*v[2]]
}