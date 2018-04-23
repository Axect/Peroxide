extern crate num_traits;

use self::num_traits::Num;

pub type Row<T> = Vec<T>;
pub type Col<T> = Vec<T>;
pub type Matrix<T> = Vec<Vec<T>>;

pub trait Add<RHS = Self> {
    type Output;

    fn add(&self, rhs: &RHS) -> Self::Output;
}

impl<T: Num + Copy> Add<T> for Matrix<T> {
    type Output = Matrix<T>;

    fn add(&self, other: &T) -> Matrix<T> {
        self.iter()
            .map(|x| x.into_iter().map(|t| *t + *other).collect::<Row<T>>())
            .collect::<Matrix<T>>()
    }
}

impl<T: Num + Copy> Add<Matrix<T>> for Matrix<T> {
    type Output = Matrix<T>;

    fn add(&self, other: &Matrix<T>) -> Matrix<T> {
        assert!(
            (self.len(), self[0].len()) == (other.len(), (other[0].len())),
            "Length does not match!"
        );
        let (m, n) = (self.len(), self[0].len());
        let mut a: Matrix<T> = vec![vec![T::zero(); n]; m];
        for i in 0..m {
            for j in 0..n {
                a[i][j] = self[i][j] + other[i][j];
            }
        }
        a
    }
}

pub trait LinAlg<T> {
    fn col(&self, index: usize) -> Col<T>;
    fn diag(&self) -> Row<T>;
    fn trace(&self) -> T;
    fn transpose(&self) -> Matrix<T>;
}

impl<T: Num + Copy> LinAlg<T> for Matrix<T> {
    fn col(&self, index: usize) -> Col<T> {
        let mut a: Col<T> = Vec::new();
        for i in 0..self.len() {
            a.push(self[i][index]);
        }
        a
    }

    fn diag(&self) -> Row<T> {
        let mut a: Row<T> = Vec::new();
        for i in 0..self.len() {
            a.push(self[i][i]);
        }
        a
    }

    fn trace(&self) -> T {
        let mut s: T = T::zero();
        for i in 0..self.len() {
            s = s + self[i][i];
        }
        s
    }

    fn transpose(&self) -> Matrix<T> {
        let m = self.len();
        let n = self[0].len();
        let mut a: Matrix<T> = vec![vec![T::zero(); m]; n];
        for i in 0..m {
            for j in 0..n {
                a[j][i] = self[i][j];
            }
        }
        a
    }
}
