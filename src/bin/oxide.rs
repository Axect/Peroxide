extern crate peroxide;

use peroxide::*;

fn main() {
    let qr = quot_rem(5, 3);
    println!("{}, {}", qr.0, qr.1);
    let e = matrix!(1;9;1, 3, 3, Row);
    println!("{}", e);
    let m = e.block();
    println!("{}", m.0);
    println!("{}", m.1);
    println!("{}", m.2);
    println!("{}", m.3);
}