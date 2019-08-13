extern crate peroxide;
use peroxide::*;

fn main() {
    let a = rand(1000, 1000);
    let pqlu = a.lu().unwrap();
}
