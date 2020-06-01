extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = Matrix::read("example_data/test.csv", false, ',');
    match a {
        Ok(m) => m.print(),
        Err(err) => println!("{}", err),
    }
}
