extern crate peroxide;

use peroxide::*;

fn main() {
    let mut a = matrix!(1;4;1, 2, 2, Row);
    a[(1,1)] = 10f64;
    println!("{}", a);
}
