mod linalg;

use linalg::Matrix;
use linalg::LinAlg;

fn main() {
    let a: Matrix<i64> = vec![vec![1,2], vec![3,4]];
    let b: Matrix<i64> = vec![vec![-1, -2], vec![-3,-4]];
    println!("{:?}", a.trace());
    println!("{:?}", a);
    println!("{:?}", a.transpose());
    println!("{:?}", a.col(1));
    println!("{:?}", a.diag());
    println!("{:?}", a.add(&b));
}
