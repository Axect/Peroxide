extern crate peroxide;
use peroxide::*;

fn main() {
    let a = zeros(1200, 1200);
    let _b = matmul(&a, &a);
}

fn matmul(a: &Matrix, b: &Matrix) -> Matrix {
    match (a.row, a.col) {
        (p, q) if p <= 100 && q <= 100 => a * b,
        _ => {
            let (a1, a2, a3, a4) = a.block();
            let (b1, b2, b3, b4) = b.block();

            let m1 = matmul(&a1,&b1) + matmul(&a2,&b3);
            let m2 = matmul(&a1,&b2) + matmul(&a2,&b4);
            let m3 = matmul(&a3,&b1) + matmul(&a4,&b3);
            let m4 = matmul(&a3,&b2) + matmul(&a4,&b4);

            combine(m1, m2, m3, m4)
        }
    }
}