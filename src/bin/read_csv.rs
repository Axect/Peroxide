extern crate peroxide;
use peroxide::*;

fn main() {
    let a = Matrix::read("example_data/test.csv", false, ',');
    match a {
        Ok(m) => m.print(),
        Err(err) => println!("{}", err),
    }
}
