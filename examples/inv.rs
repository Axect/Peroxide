extern crate peroxide;
use peroxide::*;

fn main() {
    let a = ml_matrix("1 2;3 4");
    a.det().print();
    let opt_dgrf = lapack_dgetrf(&a);
    match opt_dgrf {
        Some(dgrf) => {
            let a_inv = lapack_dgetri(&dgrf);
            match a_inv {
                Some(m) => m.print(),
                None => println!("None2"),
            }
        }
        None => {
            println!("None")
        }
    }
}