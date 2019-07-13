extern crate peroxide;
use peroxide::*;

fn main() {
    let a = Matrix::read("test.csv", false, ',');
    match a {
        Ok(m) => m.print(),
        Err(err) => println!("{}", err),
    }
}
