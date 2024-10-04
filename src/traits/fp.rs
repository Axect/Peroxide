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

/// Functional Programming for Matrix and ComplexMatrix
pub trait FPMatrix {
    type Scalar;

    fn take_row(&self, n: usize) -> Self;
    fn take_col(&self, n: usize) -> Self;
    fn skip_row(&self, n: usize) -> Self;
    fn skip_col(&self, n: usize) -> Self;
    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> Self::Scalar;
    fn col_map<F>(&self, f: F) -> Self
    where
        F: Fn(Vec<Self::Scalar>) -> Vec<Self::Scalar>;
    fn row_map<F>(&self, f: F) -> Self
    where
        F: Fn(Vec<Self::Scalar>) -> Vec<Self::Scalar>;
    fn col_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<Self::Scalar>) -> Vec<Self::Scalar>;
    fn row_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Vec<Self::Scalar>) -> Vec<Self::Scalar>;
    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
        T: Into<Self::Scalar>;
    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar;
    fn col_reduce<F>(&self, f: F) -> Vec<Self::Scalar>
    where
        F: Fn(Vec<Self::Scalar>) -> Self::Scalar;
    fn row_reduce<F>(&self, f: F) -> Vec<Self::Scalar>
    where
        F: Fn(Vec<Self::Scalar>) -> Self::Scalar;
}
