use std::ops::{Index, IndexMut};

pub trait MatrixTrait: Sized + Index<(usize, usize), Output=Self::Scalar> + IndexMut<(usize, usize)> {
    type Scalar;
    fn ptr(&self) -> *const Self::Scalar; 
    fn mut_ptr(&mut self) -> *mut Self::Scalar;
    fn as_slice(&self) -> &[Self::Scalar];
    fn as_mut_slice(&mut self) -> &mut [Self::Scalar];
    fn change_shape(&self) -> Self;
    fn change_shape_mut(&mut self);
    fn spread(&self) -> String;
    fn col(&self, index: usize) -> Vec<Self::Scalar>;
    fn row(&self, index: usize) -> Vec<Self::Scalar>;
    fn diag(&self) -> Vec<Self::Scalar>;
    fn transpose(&self) -> Self;
    fn t(&self) -> Self { self.transpose() }
    fn subs_col(&mut self, idx: usize, v: &[Self::Scalar]);
    fn subs_row(&mut self, idx: usize, v: &[Self::Scalar]);
    fn from_index<F, G>(f: F, size: (usize, usize)) -> Self
    where
        F: Fn(usize, usize) -> G + Copy,
        G: Into<Self::Scalar>;
    fn to_vec(&self) -> Vec<Vec<Self::Scalar>>;
    fn to_diag(&self) -> Self;
    fn submat(&self, start: (usize, usize), end: (usize, usize)) -> Self;
    fn subs_mat(&mut self, start: (usize, usize), end: (usize, usize), m: &Self);
}

// ┌─────────────────────────────────────────────────────────┐
//  For Linear Algebra
// └─────────────────────────────────────────────────────────┘
/// Linear algebra trait
pub trait LinearAlgebra<M:MatrixTrait> {
    fn back_subs(&self, b: &[M::Scalar]) -> Vec<M::Scalar>;
    fn forward_subs(&self, b: &[M::Scalar]) -> Vec<M::Scalar>;
    fn lu(&self) -> PQLU<M>;
    fn waz(&self, d_form: Form) -> Option<WAZD<M>>;
    fn qr(&self) -> QR<M>;
    fn svd(&self) -> SVD<M>;
    #[cfg(feature = "O3")]
    fn cholesky(&self, uplo: UPLO) -> M;
    fn rref(&self) -> M;
    fn det(&self) -> M::Scalar;
    fn block(&self) -> (M, M, M, M);
    fn inv(&self) -> M;
    fn pseudo_inv(&self) -> M;
    fn solve(&self, b: &[M::Scalar], sk: SolveKind) -> Vec<M::Scalar>;
    fn solve_mat(&self, m: &M, sk: SolveKind) -> M;
    fn is_symmetric(&self) -> bool;
}

#[allow(non_snake_case)]
pub fn solve<M:MatrixTrait + LinearAlgebra<M>>(A: &M, b: &M, sk: SolveKind) -> M {
    A.solve_mat(b, sk)
}

/// Data structure for Complete Pivoting LU decomposition
///
/// # Usage
/// ```rust
/// use peroxide::fuga::*;
///
/// let a = ml_matrix("1 2;3 4");
/// let pqlu = a.lu();
/// let (p, q, l, u) = pqlu.extract();
/// // p, q are permutations
/// // l, u are matrices
/// l.print(); // lower triangular
/// u.print(); // upper triangular
/// ```
#[derive(Debug, Clone)]
pub struct PQLU<M: MatrixTrait> {
    pub p: Vec<usize>,
    pub q: Vec<usize>,
    pub l: M,
    pub u: M,
}

#[derive(Debug, Clone)]
pub struct WAZD<M: MatrixTrait> {
    pub w: M,
    pub z: M,
    pub d: M,
}

#[derive(Debug, Clone)]
pub struct QR<M: MatrixTrait> {
    pub q: M,
    pub r: M,
}

#[derive(Debug, Copy, Clone)]
pub enum Form {
    Diagonal,
    Identity,
}

#[derive(Debug, Copy, Clone)]
pub enum SolveKind {
    LU,
    WAZ,
}

impl<M:MatrixTrait> QR<M> {
    pub fn q(&self) -> &M {
        &self.q
    }

    pub fn r(&self) -> &M {
        &self.r
    }
}

#[derive(Debug, Clone)]
pub struct SVD<M: MatrixTrait> {
    pub s: Vec<f64>,
    pub u: M,
    pub vt: M,
}

