extern crate peroxide;
use peroxide::*;
use self::Header::*;
use std::fmt;

fn main() {
    let a = DataFrame::new(vec![X, Y], vec![c!(1,2,3), c!(4,5,6)]);
    println!("{}", a);
}

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub enum Header {
    X,
    Y
}
