extern crate peroxide;

use peroxide::*;

fn main() {
    let a = matrix!(1;4;1, 2, 2, Col);
    let pqlu = a.lu();
    match pqlu {
        None => println!("Singular"),
        Some(m) => {
            let (p, q, l, u) = (m.p, m.q, m.l, m.u);
            println!("{:?}", p);
            println!("{:?}", q);
            println!("{}", l);
            println!("{}", u);
        }
    }
}
