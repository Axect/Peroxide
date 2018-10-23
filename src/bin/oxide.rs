extern crate peroxide;

use peroxide::*;

fn main() {
    let a = matrix!(1;4;1, 2, 2, Row);
    let (p,l,u) = a.plu();
    println!("{}\n{}\n{}", p, l ,u);
}