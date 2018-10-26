extern crate peroxide;

use peroxide::*;

fn main() {
    let a = matrix(c!(1,2,2,7,1,5,3,4,2), 3, 3, Row);
    let pqlu = a.lu();
    match pqlu {
        None => println!("None!"),
        Some(m) => {
            println!("{:?}\n{:?}\n{}\n{}", m.p, m.q, m.l, m.u);
            println!("{}", a.det());
            println!("{}", a.inv().unwrap());
        }
    }
}