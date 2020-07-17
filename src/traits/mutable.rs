use crate::structure::matrix::Shape;

pub trait MutFP {
    type Scalar;
    fn mut_map<F>(&mut self, f: F)
    where
        F: Fn(Self::Scalar) -> Self::Scalar;
    fn mut_zip_with<F>(&mut self, f: F, other: &Self)
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar;
}

pub trait MutMatrix {
    unsafe fn col_mut(&mut self, idx: usize) -> Vec<*mut f64>;
    unsafe fn row_mut(&mut self, idx: usize) -> Vec<*mut f64>;
    unsafe fn swap(&mut self, idx1: usize, idx2: usize, shape: Shape);
    unsafe fn swap_with_perm(&mut self, p: &Vec<(usize, usize)>, shape: Shape);
}
