pub type Row = Vec<u64>;
pub type Col = Vec<u64>;
pub type Matrix = Vec<Vec<u64>>;

pub trait LinAlg {
    fn trace(&self) -> u64;
    fn transpose(&self) -> Matrix;
}

impl LinAlg for Matrix {
    fn trace(&self) -> u64 {
        let mut s: u64 = 0;
        for i in 0..self.len() {
            s += self[i][i];
        }
        s
    }

    fn transpose(&self) -> Matrix {
        let m = self.len();
        let n = self[0].len();
        let mut a: Matrix = vec![vec![0; m]; n];
        for i in 0..m {
            for j in 0..n {
                a[j][i] = self[i][j];
            }
        }
        a
    }
}
