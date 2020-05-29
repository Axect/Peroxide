use std::convert::Into;

/// Mathematical Vector
///
/// # Description
/// Vector has two operations : addition, scalar multiplication.
/// And a space of the vector should closed for that operations.
pub trait Vector {
    fn add_vec<'a, 'b>(&'a self, rhs: &'b Self) -> Self;
    fn sub_vec<'a, 'b>(&'a self, rhs: &'b Self) -> Self;
    fn mul_scalar<T: Into<f64>>(&self, rhs: T) -> Self;
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
    type Scalar;
    fn norm(&self, kind: Norm) -> Self::Scalar;
}

/// Inner product Vector
pub trait InnerProduct: Normed {
    fn dot(&self, rhs: &Self) -> Self::Scalar;
}

/// Linear operation for Vector
pub trait LinearOp<T: Vector> {
    fn apply(&self, rhs: &T) -> T;
}