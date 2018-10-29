extern crate peroxide;

use peroxide::*;

fn main() {
    let a = matrix!(1;6;1, 3, 2, Row);
    let blocks = a.block();
    println!("{}", blocks.0);
    println!("{}", blocks.1);
    println!("{}", blocks.2);
    println!("{}", blocks.3);
}
