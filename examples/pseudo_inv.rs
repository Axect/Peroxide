extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = ml_matrix("2 -1 0;4 3 -2");
    #[cfg(feature="O3")]
    {
        let b = a.pseudo_inv();
        assert_eq!(&a * &b, eye(2));
    }
    a.print();
}
