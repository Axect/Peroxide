extern crate peroxide;
use peroxide::*;

fn main() {
    let a = ml_matrix("1 2;3 4");
    let (q, r) = modified_gram_schmidt(a.clone());
    q.print();
    r.print();
    (q % r).print();
}

fn modified_gram_schmidt(mat: Matrix) -> (Matrix, Matrix) {
    let m = mat.row;
    let n = mat.col;

    let mut a = mat.clone();
    let mut q = zeros(m,n);
    let mut r = zeros(n,n);

    for k in 0 .. n {
        let mut s = 0f64;
        for j in 0 .. m {
            s += a[(j, k)].powi(2);
        }
        r[(k, k)] = s.sqrt();
        for j in 0 .. m {
            q[(j, k)] = a[(j, k)] / r[(k, k)];
        }
        for i in (k+1) .. n {
            s = 0f64;
            for j in 0 .. m {
                s += a[(j,i)] * a[(j,k)];
            }
            r[(k, i)] = s;
            for j in 0 .. m {
                a[(j,i)] -= r[(k,i)] * q[(j,k)];
            }
        }
    }
    (q, r)
}