extern crate peroxide;

use peroxide::*;

fn main() {
    let a = matrix(c!(1,2,2,4,5,1,7,2,9), 3, 3, Row);
    let (l, u) = a.lu();
    println!("{}\n{}", l, u);

    let mut m = a.data.clone();
    m[0] -= 1.;
    m[4] += 3.;
    m[8] -= 23.;
    let n = matrix(m, 3, 3, Row);
    println!("{}", n);
    let (l2, u2) = n.lu();
    println!("{}\n{}", l2,u2);
}