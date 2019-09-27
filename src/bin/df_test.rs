extern crate peroxide;
use peroxide::*;
use self::Header::*;
use std::fmt;

fn main() {
    let a = DataFrame::new(vec![X, Y, Z], vec![c!(1,2,3), c!(4,5,6), c!(7,8,9)]);
    println!("{}", a);
}

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub enum Header {
    X,
    Y,
    Z
}
