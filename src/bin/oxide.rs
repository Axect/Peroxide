extern crate peroxide;

use peroxide::*;

fn main() {
    let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
    println!("{}", a);
    let b = Matrix::new(vec![1,2,3,4], 2, 2, Col);
    println!("{}", b);
    println!("{}", a.clone() * b.clone());
    println!("{}", a.clone() % b.clone());

    let c = seq!(1,10,1);
    println!("{:?}", c);

    println!("{}", a);
    println!("{:?}", a.col(0));
    println!("{:?}", a.col(1));
    println!("{:?}", a.row(0));
    println!("{:?}", a.row(1));
    println!("{}", a[(0,0)]);
    println!("{}", a[(0,1)]);
    println!("{}", a[(1,0)]);
    println!("{}", a[(1,1)]);

    let lu = a.lu();
    println!("{}", lu.0);
    println!("{}", lu.1);

    let d = matrix(vec![1,2,2,4,1,5,6,-2,3], 3, 3, Row);
    let lu2 = d.lu();
    println!("{}", lu2.0);
    println!("{}", lu2.1);
    println!("{}", lu2.0 % lu2.1);

    let e = matrix!(1;16;1, 4, 4, Row);
    let m = e.block();
    println!("{}", m.0);
    println!("{}", m.1);
    println!("{}", m.2);
    println!("{}", m.3);
}