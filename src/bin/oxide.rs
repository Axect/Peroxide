extern crate peroxide;

use peroxide::*;

fn main() {
    let a = matrix!(1;10000;1, 100, 100, Row);
    println!("{}", a.clone() % a);
}
