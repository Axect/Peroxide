mod linalg;

use linalg::Matrix;
use linalg::LinAlg;

fn main() {
    let a: Matrix<i64> = vec![vec![1,2], vec![3,4]];
    println!("{:?}", a.trace());
    println!("{:?}", a);
    println!("{:?}", a.transpose());
}
