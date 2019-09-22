extern crate peroxide;
use peroxide::*;

fn main() {
    let mut a = zeros(10, 100_0000);
    a.swap(0, 9, Col);
    a.print();
}
