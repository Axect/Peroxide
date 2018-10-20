extern crate peroxide;

use peroxide::*;

fn main() {
    let a = Matrix::new(vec![1,2,3,4], 2, 2, Row);
    println!("{}", a);
    let b = Matrix::new(vec![1,2,3,4], 2, 2, Col);
    println!("{}", b);
    // println!("{}", a.clone() + b.clone());
    // println!("{}", a.clone() - b.clone());
    println!("{}", a.clone() * b.clone());
    // println!("{}", a.fmap(|x| x + 1f64));
    // println!("{}", a.reduce(1, |x,y| x*y));
    println!("{}", a.zip_with(|x,y| x * y, &b));
    // println!("{}", a.clone() + 1f64);
    println!("{}", a.clone() % b.clone());

    let c = seq(1,10,1);
    println!("{}", c);

    println!("{}", a);
    println!("{}", a.col(0));
    println!("{}", a.col(1));
    println!("{}", a.row(0));
    println!("{}", a.row(1));
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

    let f = seq(1, 10, 1);
    println!("{}", f);
}