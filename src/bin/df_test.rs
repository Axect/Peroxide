extern crate peroxide;
use peroxide::*;
use self::Header::*;
use std::fmt;

fn main() {
    let mut a: DataFrame<Header> = DataFrame::new();
    a.insert(X, c!(1,2,3));
    a.insert(Y, c!(4,5,6));
    a.insert(Z, c!(7,8,9));
    println!("{}", a);

    let b: DataFrame<Header> = DataFrame::from_matrix(vec![X, Y, Z], a.to_matrix());
    println!("{}", b);
}

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub enum Header {
    X,
    Y,
    Z
}
