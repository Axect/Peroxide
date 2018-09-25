extern crate peroxide;

use peroxide::*;

fn main() {
    let mut a = Matrix::new(vec![1,2,3,4], 2, 2, Col);
    println!("{}", a);
    a = a.change_shape();
    println!("{}", a.shape);
    println!("{}", a.transpose());
}