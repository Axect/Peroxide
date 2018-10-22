extern crate peroxide;

use peroxide::*;

fn main() {
    let e = matrix!(1;9;1, 3, 3, Row);
    println!("{}", e);
    let m = e.block();
    println!("{}", m.0);
    println!("{}", m.1);
    println!("{}", m.2);
    println!("{}", m.3);

    let f = matrix!(1;9;1, 3, 3, Col);
    println!("{}", f);
    let n = f.block();
    println!("{}", n.0);
    println!("{}", n.1);
    println!("{}", n.2);
    println!("{}", n.3);

    let mc = combine(m.0, m.1, m.2, m.3);
    let nc = combine(n.0, n.1, n.2, n.3);

    println!("{}", mc);
    println!("{}", nc);
}