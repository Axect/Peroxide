#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;
use peroxide::structure::dataframe::*;
use peroxide::structure::dataframe::DTypeArray::*;

fn main() {
    let a = Series::new(USIZE(vec![1,2,3]));
    let b = Series::new(Str(vec!["i", "j", "k"].into_iter().map(|x| x.into()).collect()));

    let mut df = DataFrame::new(vec![a, b]);

    println!("{:?}", df["0"]);
    println!("{:?}", df["1"]);

    df["1"] = Series::new(F64(c!(1,2,3)));

    println!("{:?}", df[1]);

    df.push("a", Series::new(Char(vec!['a', 'b', 'c'])));

    println!("{:?}", df);
    println!("{:?}", df.row(1));

    println!("{}", df.row(1)["a"].at(0).unwrap_char());
}
