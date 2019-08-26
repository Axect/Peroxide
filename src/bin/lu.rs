extern crate peroxide;
use peroxide::*;

fn main() {
    #[cfg(feature = "openblas")]
    {
        let a = ml_matrix("1 2;3 4");
        let opt_dgrf = lapack_dgetrf(&a);
        match opt_dgrf {
            None => println!("None"),
            Some(dgrf) => {
                let l = dgrf.get_L();
                let u = dgrf.get_U();
                l.print();
                u.print();
            }
        }
    }
}