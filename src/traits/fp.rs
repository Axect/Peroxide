use std::convert::Into;
use crate::structure::matrix::Matrix;

/// Functional Programming tools for Vector
pub trait FPVector {
    type Scalar;

    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> Self::Scalar;
    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
        T: Into<Self::Scalar>;
    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar;
    fn filter<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> bool;
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
        F: Fn(f64) -> f64;
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
        F: Fn(f64, f64) -> f64,
        T: Into<f64>;
    fn zip_with<F>(&self, f: F, other: &Matrix) -> Matrix
    where
        F: Fn(f64, f64) -> f64;
}
