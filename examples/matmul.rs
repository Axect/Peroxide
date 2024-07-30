extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let row = 1000usize;
    let col = 1000usize;

    // Create Matrix
    let m = rand(row, col);

    // Create another Matrix
    let n = rand(row, col);

    // Matmul
    let result = m * n;

    result[(row / 2, col / 2)].print();
}
