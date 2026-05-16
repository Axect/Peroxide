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
    type Scalar;

    /// Returns raw mutable pointers to the elements of column `idx`.
    ///
    /// # Safety
    ///
    /// `idx` must be a valid column index for the matrix and the
    /// caller must ensure that no other reference to the matrix is
    /// used concurrently while the returned pointers are alive. The
    /// pointers are invalidated if the matrix is reallocated.
    unsafe fn col_mut(&mut self, idx: usize) -> Vec<*mut Self::Scalar>;

    /// Returns raw mutable pointers to the elements of row `idx`.
    ///
    /// # Safety
    ///
    /// `idx` must be a valid row index for the matrix and the caller
    /// must ensure that no other reference to the matrix is used
    /// concurrently while the returned pointers are alive. The
    /// pointers are invalidated if the matrix is reallocated.
    unsafe fn row_mut(&mut self, idx: usize) -> Vec<*mut Self::Scalar>;

    /// Swap row or column `idx1` with `idx2` in place.
    ///
    /// # Safety
    ///
    /// `idx1` and `idx2` must both be valid indices for the given
    /// `shape` (rows or columns) of the matrix.
    unsafe fn swap(&mut self, idx1: usize, idx2: usize, shape: Shape);

    /// Apply a permutation of (idx1, idx2) swaps in order.
    ///
    /// # Safety
    ///
    /// Every `(i, j)` in `p` must be a valid index pair for the given
    /// `shape` of the matrix.
    unsafe fn swap_with_perm(&mut self, p: &[(usize, usize)], shape: Shape);
}

// Mutable trait for Vector in Parallel (Uses Rayon crate)
pub trait ParallelMutFP {
    type Scalar;
    fn par_mut_map<F>(&mut self, f: F)
    where
        F: Fn(Self::Scalar) -> Self::Scalar + Send + Sync;
    fn par_mut_zip_with<F>(&mut self, f: F, other: &Self)
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar + Send + Sync;
}
