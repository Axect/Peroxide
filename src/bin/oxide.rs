extern crate peroxide;

use peroxide::*;

fn main() {
    let a = Matrix::new(vec![1,2,3,4], 2, 2, Col);
    println!("{}", a);
}