//! Sparse matrix (CCS format)
//!
//! * Reference : Press, William H., and William T. Vetterling. *Numerical Recipes.* Cambridge: Cambridge Univ. Press, 2007.

use crate::structure::matrix::{Form, LinearAlgebra, Matrix, SolveKind, PQLU, QR, WAZD};
use crate::traits::math::{InnerProduct, LinearOp, Norm, Normed, Vector};
use crate::util::non_macro::zeros;
use std::ops::Mul;

#[derive(Debug, Clone)]
pub struct SPMatrix {
    pub row: usize,
    pub col: usize,
    pub nnz: usize,
    pub col_ptr: Vec<usize>,
    pub row_ics: Vec<usize>,
    pub data: Vec<f64>,
}

impl SPMatrix {
    pub fn new(row: usize, col: usize, nnz: usize) -> Self {
        SPMatrix {
            row,
            col,
            nnz,
            col_ptr: vec![0usize; col + 1],
            row_ics: vec![0usize, nnz],
            data: vec![0f64; nnz],
        }
    }

    pub fn from_dense(m: &Matrix) -> Self {
        let mut data: Vec<f64> = Vec::new();
        let mut row_ics: Vec<usize> = Vec::new();
        let mut col_ptr: Vec<usize> = vec![0usize; m.col + 1];
        let mut k = 0usize;
        for j in 0..m.col {
            for i in 0..m.row {
                let val = m[(i, j)];
                if val != 0f64 {
                    data.push(val);
                    row_ics.push(i);
                    k += 1;
                }
            }
            col_ptr[j + 1] = k;
        }

        SPMatrix {
            row: m.row,
            col: m.col,
            nnz: data.len(),
            col_ptr,
            row_ics,
            data,
        }
    }

    pub fn to_dense(&self) -> Matrix {
        let mut m = zeros(self.row, self.col);

        for j in 0..self.col {
            for i in self.col_ptr[j]..self.col_ptr[j + 1] {
                let k = self.row_ics[i];
                m[(k, j)] = self.data[i];
            }
        }
        m
    }

    pub fn col_ptr(&self) -> &Vec<usize> {
        &self.col_ptr
    }

    pub fn row_ics(&self) -> &Vec<usize> {
        &self.row_ics
    }

    pub fn data(&self) -> &Vec<f64> {
        &self.data
    }

    pub fn transpose(&self) -> Self {
        let row = self.row;
        let col = self.col;
        let nnz = self.nnz;
        let col_ptr = self.col_ptr();
        let row_ics = self.row_ics();
        let data = self.data();
        let mut count = vec![0usize; row];
        let mut result = Self::new(col, row, nnz);

        for i in 0..col {
            for j in col_ptr[i]..col_ptr[i + 1] {
                let k = row_ics[j];
                count[k] += 1;
            }
        }
        for j in 0..row {
            result.col_ptr[j + 1] = result.col_ptr[j] + count[j];
            count[j] = 0;
        }
        for i in 0..col {
            for j in col_ptr[i]..col_ptr[i + 1] {
                let k = row_ics[j];
                let index = result.col_ptr[k] + count[k];
                result.row_ics[index] = i;
                result.data[index] = data[j];
                count[k] += 1;
            }
        }
        result
    }

    pub fn t(&self) -> Self {
        self.transpose()
    }
}

impl LinearOp<Vec<f64>, Vec<f64>> for SPMatrix {
    fn apply(&self, rhs: &Vec<f64>) -> Vec<f64> {
        let mut y = vec![0f64; self.row];
        let col_ptr = self.col_ptr();
        let row_ics = self.row_ics();
        let data = self.data();
        for j in 0..self.col {
            for i in col_ptr[j]..col_ptr[j + 1] {
                y[row_ics[i]] += data[i] * rhs[j];
            }
        }
        y
    }
}

/// Linear algebra for sparse matrix
///
/// **Caution** : In every ops in this trait, there is converting process to dense matrix
impl LinearAlgebra for SPMatrix {
    fn back_subs(&self, _b: &Vec<f64>) -> Vec<f64> {
        unimplemented!()
    }

    fn forward_subs(&self, _b: &Vec<f64>) -> Vec<f64> {
        unimplemented!()
    }

    fn lu(&self) -> PQLU {
        self.to_dense().lu()
    }

    fn waz(&self, _d_form: Form) -> Option<WAZD> {
        unimplemented!()
    }

    fn qr(&self) -> QR {
        self.to_dense().qr()
    }

    fn det(&self) -> f64 {
        self.to_dense().det()
    }

    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix) {
        self.to_dense().block()
    }

    fn inv(&self) -> Matrix {
        self.to_dense().inv()
    }

    fn pseudo_inv(&self) -> Matrix {
        self.to_dense().pseudo_inv()
    }

    fn rref(&self) -> Matrix {
        self.to_dense().rref()
    }

    fn solve(&self, _b: &Vec<f64>, _sk: SolveKind) -> Vec<f64> {
        unimplemented!()
    }

    fn solve_mat(&self, _m: &Matrix, _sk: SolveKind) -> Matrix {
        unimplemented!()
    }
}

/// Matrix multiplication with vector
impl Mul<Vec<f64>> for SPMatrix {
    type Output = Vec<f64>;
    fn mul(self, rhs: Vec<f64>) -> Self::Output {
        self.apply(&rhs)
    }
}

/// Reference version of matrix multiplication with vector
impl<'a, 'b> Mul<&'b Vec<f64>> for &'a SPMatrix {
    type Output = Vec<f64>;
    fn mul(self, rhs: &'b Vec<f64>) -> Self::Output {
        self.apply(rhs)
    }
}

impl Into<Matrix> for SPMatrix {
    fn into(self) -> Matrix {
        self.to_dense()
    }
}

impl Into<SPMatrix> for Matrix {
    fn into(self) -> SPMatrix {
        SPMatrix::from_dense(&self)
    }
}
