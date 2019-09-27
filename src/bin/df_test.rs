extern crate peroxide;
use peroxide::*;
use self::Header::*;

fn main() {
    let mut a: DataFrame<Header> = DataFrame::new();
    a.insert(X, seq(0, 100, 1));
    a.insert(Y, c!(4,5,6));
    a.insert(Z, c!(7,8,9));
    a.print();
}

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub enum Header {
    X,
    Y,
    Z
}
