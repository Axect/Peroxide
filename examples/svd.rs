extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = ml_matrix("14 2;4 22;16 13") / 15f64;
    #[cfg(feature="O3")]
    {
        let svd = a.svd();
        svd.u.print();
        svd.vt.print();
        svd.s.print();
    }
    a.print();
}
