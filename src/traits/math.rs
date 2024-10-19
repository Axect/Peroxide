use crate::structure::matrix::Matrix;

/// Mathematical Vector
///
/// # Description
/// Vector has two operations : addition, scalar multiplication.
/// And a space of the vector should closed for that operations.
pub trait Vector {
    type Scalar;
    fn add_vec(&self, rhs: &Self) -> Self;
    fn sub_vec(&self, rhs: &Self) -> Self;
    fn mul_scalar(&self, rhs: Self::Scalar) -> Self;
}

/// Kinds of Vector & Matrix norm
///
/// # Kinds of Vector norm
/// * `l1`
/// * `l2`
/// * `lp`
/// * `lInf`
///
/// # Kinds of Matrix norm
/// * `F`: Frobenius norm
/// * `lpq`: Element-wise pq norm
#[derive(Debug, Copy, Clone)]
pub enum Norm {
    L1,
    L2,
    Lp(f64),
    LInf,
    F,
    Lpq(f64, f64),
}

/// Normed Vector
pub trait Normed: Vector {
    type UnsignedScalar;
    fn norm(&self, kind: Norm) -> Self::UnsignedScalar;
    fn normalize(&self, kind: Norm) -> Self
    where
        Self: Sized;
}

/// Inner product Vector
pub trait InnerProduct: Normed {
    fn dot(&self, rhs: &Self) -> Self::Scalar;
}

/// Linear operation for Vector
pub trait LinearOp<T: Vector, S: Vector> {
    fn apply(&self, rhs: &T) -> S;
}

/// Vector Products
pub trait VectorProduct: Vector {
    fn cross(&self, other: &Self) -> Self;
    fn outer(&self, other: &Self) -> Matrix;
}

/// Matrix Products
pub trait MatrixProduct {
    fn kronecker(&self, other: &Self) -> Self;
    fn hadamard(&self, other: &Self) -> Self;
}

// =============================================================================
// Implementation for primitive types
// =============================================================================

impl Vector for f64 {
    type Scalar = Self;

    fn add_vec(&self, rhs: &Self) -> Self {
        self + rhs
    }

    fn sub_vec(&self, rhs: &Self) -> Self {
        self - rhs
    }

    fn mul_scalar(&self, rhs: Self::Scalar) -> Self {
        self * rhs
    }
}

impl Normed for f64 {
    type UnsignedScalar = f64;
    fn norm(&self, _kind: Norm) -> Self::Scalar {
        self.abs()
    }

    fn normalize(&self, _kind: Norm) -> Self
    where
        Self: Sized,
    {
        self / self.abs()
    }
}

// =============================================================================
// Implementation for parallel traits
// =============================================================================

/// Mathematical Vector in Parallel
pub trait ParallelVector {
    type Scalar;
    fn par_add_vec(&self, rhs: &Self) -> Self;
    fn par_sub_vec(&self, rhs: &Self) -> Self;
    fn par_mul_scalar(&self, rhs: Self::Scalar) -> Self;
}

/// Normed Vector in Parallel
pub trait ParallelNormed: Vector {
    type UnsignedScalar;
    fn par_norm(&self, kind: Norm) -> Self::UnsignedScalar;
}

/// Inner product Vector in Parallel
pub trait ParallelInnerProduct: ParallelNormed {
    fn par_dot(&self, rhs: &Self) -> Self::Scalar;
}

/// Matrix Products in Parallel
pub trait ParallelMatrixProduct {
    fn par_hadamard(&self, other: &Self) -> Matrix;
}

/// Vector Products in Parallel
pub trait ParallelVectorProduct: Vector {
    fn par_cross(&self, other: &Self) -> Self;
}
