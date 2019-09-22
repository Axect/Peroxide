extern crate peroxide;
use peroxide::*;

fn main() {
    let mut a = zeros(10, 100000);
    unsafe {
        let mut p = a.col_mut(0);
        for ptr in p {
            *ptr = 1f64;
        }
        a.col(0)[0].print();
    }
}