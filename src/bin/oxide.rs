extern crate peroxide;

use peroxide::*;

fn main() {
    let m = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    println!("{:?}", m.var());
    println!("{:?}", m.sd());
}
