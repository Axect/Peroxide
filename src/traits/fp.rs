use crate::structure::matrix::Matrix;

/// Functional Programming tools for Vector
pub trait FPVector {
    type Scalar;

    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> Self::Scalar + Send + Sync;
    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar + Send + Sync,
        T: Into<Self::Scalar> + Send + Sync + Copy;
    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar + Send + Sync;
    fn filter<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> bool + Send + Sync;
    fn take(&self, n: usize) -> Self;
    fn skip(&self, n: usize) -> Self;
    fn sum(&self) -> Self::Scalar;
    fn prod(&self) -> Self::Scalar;
}

/// Functional Programming for Matrix
pub trait FPMatrix {
    fn take_row(&self, n: usize) -> Matrix;
    fn take_col(&self, n: usize) -> Matrix;
    fn skip_row(&self, n: usize) -> Matrix;
    fn skip_col(&self, n: usize) -> Matrix;
    fn fmap<F>(&self, f: F) -> Matrix
    where
        F: Fn(f64) -> f64 + Send + Sync;
    fn col_map<F>(&self, f: F) -> Matrix
    where
        F: Fn(Vec<f64>) -> Vec<f64>;
    fn row_map<F>(&self, f: F) -> Matrix
    where
        F: Fn(Vec<f64>) -> Vec<f64>;
    fn col_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<f64>) -> Vec<f64>;
    fn row_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<f64>) -> Vec<f64>;
    fn reduce<F, T>(&self, init: T, f: F) -> f64
    where
        F: Fn(f64, f64) -> f64 + Send + Sync,
        T: Into<f64> + Copy + Clone + Send + Sync;
    fn zip_with<F>(&self, f: F, other: &Matrix) -> Matrix
    where
        F: Fn(f64, f64) -> f64 + Send + Sync;
    fn col_reduce<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(Vec<f64>) -> f64 + Send + Sync;
    fn row_reduce<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(Vec<f64>) -> f64 + Send + Sync;
}
