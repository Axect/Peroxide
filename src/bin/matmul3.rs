extern crate peroxide;
use peroxide::*;

fn main() {
    let a = zeros(2048, 2048);
    let _b = strassen(&a, &a);
}

/// Strassen's algorithm
/// # Assumption for convenience
/// * All matrices are square matrix
/// * Except base case, size of all matrices are power of 2
fn strassen(a: &Matrix, b: &Matrix) -> Matrix {
    // Guarantee same length
    assert_eq!(a.col, b.row);

    // Base case
    match (a.row, a.col) {
        (p, q) if p <= 100 && q <= 100 => a * b,
        _ => {
            let (a1, a2, a3, a4) = a.block(); // Split a
            let (b1, b2, b3, b4) = b.block(); // Split b

            // Do Strassen's algorithm
            let p1 = strassen(&a1, &(&b2 - &b4));          // A (F - H)
            let p2 = strassen(&(&a1 + &a2), &b4);          // (A + B) H
            let p3 = strassen(&(&a3 + &a4), &b1);          // (C + D) E
            let p4 = strassen(&a4, &(&b3 - &b1));          // D (G - E)
            let p5 = strassen(&(&a1 + &a4), &(&b1 + &b4)); // (A + D)(E + H)
            let p6 = strassen(&(&a2 - &a4), &(&b3 + &b4)); // (B - D)(G + H)
            let p7 = strassen(&(&a1 - &a3), &(&b1 + &b2)); // (A - C)(E + F)

            let part1 = &(&p5 + &p4) + &(&p6 - &p2); // P5 + P4 - (P2 - P6)
            let part2 = &p1 + &p2;                       // P1 + P2
            let part3 = &p3 + &p4;                       // P3 + P4
            let part4 = &(&p1 - &p3) + &(&p5 - &p7); // P1 + P5 - (P3 + P7)

            combine(part1, part2, part3, part4)
        }
    }
}