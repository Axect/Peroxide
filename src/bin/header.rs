extern crate peroxide;
use peroxide::*;

fn main() {
    let a = m_matrix("1 2 3; 4 5 6");
    let h = vec!["a", "b", "c"];
    a.write_with_header("test_header.csv", h);
}