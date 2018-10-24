extern crate peroxide;

use peroxide::*;

fn main() {
    let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    println!("{}", a);
    println!("{:?}", a.col(0));
    println!("{:?}", a.col(1));
    println!("{}", a.cov());
    println!("{}", cov(c!(1,2,3),c!(3,2,1)));
}