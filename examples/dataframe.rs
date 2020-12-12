#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;
use peroxide::structure::dataframe::*;
use peroxide::structure::dataframe::DTypeArray::*;

fn main() {
    let x = c!(1,2,3,4);
    let a = Series::new(x);
    println!("{:?}", a);

    let s = a.at(0);
    println!("{:?}", s);

    let b = Series::new(vec!['a', 'b', 'c', 'd']);

    let mut df = DataFrame::new(vec![a, b]);

    println!("{:?}", df["0"]);
    println!("{:?}", df["1"]);

    df["1"] = Series::new(c!(5,6,7,8));

    println!("{:?}", df[1]);

    df.push("a", Series::new(vec!['a', 'b', 'c', 'd']));

    println!("{:?}", df);
    println!("{:?}", df.row(1));

    let ch: char = df.row(1)["a"].at(0).unwrap();
    println!("{}", ch);

    println!("{}", df);
}
