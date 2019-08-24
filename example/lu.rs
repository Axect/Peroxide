extern crate peroxide;
use peroxide::*;

fn main() {
    let a = ml_matrix("1 2;3 4");
    let mid = lapack_dgetrf(&a);
    match mid {
        Some(dgrf) => dgrf.fact_mat.print(),
        None => println!("None")
    }
}
