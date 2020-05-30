use crate::traits::{
    math::{Normed, Norm}
};
use crate::structure::matrix::Matrix;
use crate::numerical::{
    integral,
    integral::Integral::GaussLegendre,
    eigen,
    eigen::{Eigen, EigenMethod::Jacobi},
};

/// Simple Norm
pub trait SimpleNorm: Normed {
    fn norm(&self) -> Self::Scalar;
    fn normalize(&self) -> Self::Scalar;
}

/// Simple integrate
pub fn integrate<F: Fn(f64) -> f64>(f: F, (a, b): (f64, f64)) -> f64 {
    integral::integrate(f, (a, b), GaussLegendre(15))
}

/// Simple Eigenpair
pub fn eigen(m: &Matrix) -> Eigen {
    eigen::eigen(m, Jacobi)
}

/// Simple L2 norm
impl SimpleNorm for Vec<f64> {
    fn norm(&self) -> Self::Scalar {
        self.norm(Norm::L2)
    }

    fn normalize(&self) -> Self::Scalar {
        self.normalize(Norm::L2)
    }
}

/// Simple Frobenius norm
impl SimpleNorm for Matrix {
    fn norm(&self) -> Self::Scalar {
        self.norm(Norm::F)
    }

    fn normalize(&self) -> Self::Scalar {
        unimplemented!()
    }
}