extern crate num_traits;

use self::num_traits::Num;

pub type Row<T> = Vec<T>;
pub type Col<T> = Vec<T>;
pub type Matrix<T> = Vec<Vec<T>>;

pub trait LinAlg<T> {
    fn trace(&self) -> T;
    fn transpose(&self) -> Matrix<T>;
}

impl<T: Num + Clone> LinAlg<T> for Matrix<T> {
    fn trace(&self) -> T {
        let mut s: T = T::zero();
        for i in 0..self.len() {
            s = s + self[i][i].clone();
        }
        s
    }

    fn transpose(&self) -> Matrix<T> {
        let m = self.len();
        let n = self[0].len();
        let mut a: Matrix<T> = vec![vec![T::zero(); m]; n];
        for i in 0..m {
            for j in 0..n {
                a[j][i] = self[i][j].clone();
            }
        }
        a
    }
}
