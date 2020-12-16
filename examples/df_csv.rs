#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;
use peroxide::structure::dataframe::*;

fn main() {
    let a = Series::new(vec![1,2,3,4]);
    let b = Series::new(c!(0.1, 0.2, 0.3, 0.4));
    let c = Series::new(vec![true, false, true, false]);
    let d = Series::new(vec!['a', 'b', 'c', 'd']);

    let mut df = DataFrame::new(vec![a, b, c, d]);
    df.set_header(vec!["a", "b", "c", "d"]);

    df.print();

    df.write_csv("example_data/df_csv.csv").expect("Can't write csv");
}
