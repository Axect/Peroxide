#[cfg(feature="parquet")]
use std::error::Error;
use crate::numerical::{
    eigen,
    eigen::{Eigen, EigenMethod::Jacobi},
    integral,
    integral::Integral::G7K15R,
    spline,
    spline::{CubicHermiteSpline, SlopeMethod::Quadratic},
};
use crate::structure::matrix::{self, Matrix};
use crate::structure::polynomial;
use crate::traits::math::{Norm, Normed};
#[cfg(feature="parquet")]
use crate::structure::dataframe::{DataFrame, WithParquet};
#[cfg(feature="parquet")]
use arrow2::io::parquet::write::CompressionOptions;

/// Simple Norm
pub trait SimpleNorm: Normed {
    fn norm(&self) -> Self::Scalar;
    fn normalize(&self) -> Self;
}

/// Simple integrate
pub fn integrate<F: Fn(f64) -> f64 + Copy>(f: F, (a, b): (f64, f64)) -> f64 {
    integral::integrate(f, (a, b), G7K15R(1e-4, 20))
}

/// Simple Linear algebra
pub trait SimplerLinearAlgebra {
    fn back_subs(&self, b: &Vec<f64>) -> Vec<f64>;
    fn forward_subs(&self, b: &Vec<f64>) -> Vec<f64>;
    fn lu(&self) -> matrix::PQLU;
    fn waz_diag(&self) -> Option<matrix::WAZD>;
    fn waz(&self) -> Option<matrix::WAZD>;
    fn qr(&self) -> matrix::QR;
    #[cfg(feature="O3")]
    fn cholesky(&self) -> Matrix;
    fn rref(&self) -> Matrix;
    fn det(&self) -> f64;
    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix);
    fn inv(&self) -> Matrix;
    fn pseudo_inv(&self) -> Matrix;
    fn solve(&self, b: &Vec<f64>) -> Vec<f64>;
    fn solve_mat(&self, m: &Matrix) -> Matrix;
    fn is_symmetric(&self) -> bool;
}

/// Simple Eigenpair
pub fn eigen(m: &Matrix) -> Eigen {
    eigen::eigen(m, Jacobi)
}

/// Simple L2 norm
impl SimpleNorm for Vec<f64> {
    fn norm(&self) -> Self::Scalar {
        Normed::norm(self, Norm::L2)
    }

    fn normalize(&self) -> Self {
        Normed::normalize(self, Norm::L2)
    }
}

/// Simple Frobenius norm
impl SimpleNorm for Matrix {
    fn norm(&self) -> Self::Scalar {
        Normed::norm(self, Norm::F)
    }

    fn normalize(&self) -> Self {
        unimplemented!()
    }
}

impl SimplerLinearAlgebra for Matrix {
    fn back_subs(&self, b: &Vec<f64>) -> Vec<f64> {
        matrix::LinearAlgebra::back_subs(self, b)
    }

    fn forward_subs(&self, b: &Vec<f64>) -> Vec<f64> {
        matrix::LinearAlgebra::forward_subs(self, b)
    }

    fn lu(&self) -> matrix::PQLU {
        matrix::LinearAlgebra::lu(self)
    }

    fn waz_diag(&self) -> Option<matrix::WAZD> {
        matrix::LinearAlgebra::waz(self, matrix::Form::Diagonal)
    }

    fn waz(&self) -> Option<matrix::WAZD> {
        matrix::LinearAlgebra::waz(self, matrix::Form::Identity)
    }

    fn qr(&self) -> matrix::QR {
        matrix::LinearAlgebra::qr(self)
    }

    #[cfg(feature="O3")]
    fn cholesky(&self) -> Matrix {
        matrix::LinearAlgebra::cholesky(self, matrix::UPLO::Lower)
    }

    fn rref(&self) -> Matrix {
        matrix::LinearAlgebra::rref(self)
    }

    fn det(&self) -> f64 {
        matrix::LinearAlgebra::det(self)
    }

    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix) {
        matrix::LinearAlgebra::block(self)
    }

    fn inv(&self) -> Matrix {
        matrix::LinearAlgebra::inv(self)
    }

    fn pseudo_inv(&self) -> Matrix {
        matrix::LinearAlgebra::pseudo_inv(self)
    }

    fn solve(&self, b: &Vec<f64>) -> Vec<f64> {
        matrix::LinearAlgebra::solve(self, b, matrix::SolveKind::LU)
    }

    fn solve_mat(&self, m: &Matrix) -> Matrix {
        matrix::LinearAlgebra::solve_mat(self, m, matrix::SolveKind::LU)
    }

    fn is_symmetric(&self) -> bool {
        matrix::LinearAlgebra::is_symmetric(self)
    }
}

/// Simple solve
#[allow(non_snake_case)]
pub fn solve(A: &Matrix, m: &Matrix) -> Matrix {
    matrix::solve(A, m, matrix::SolveKind::LU)
}

/// Simple Chebyshev Polynomial (First Kind)
pub fn chebyshev_polynomial(n: usize) -> polynomial::Polynomial {
    polynomial::chebyshev_polynomial(n, polynomial::SpecialKind::First)
}

pub fn cubic_hermite_spline(node_x: &[f64], node_y: &[f64]) -> anyhow::Result<CubicHermiteSpline> {
    spline::cubic_hermite_spline(node_x, node_y, Quadratic)
}

/// Simple handle parquet
#[cfg(feature="parquet")]
pub trait SimpleParquet: Sized {
    fn write_parquet(&self, path: &str) -> Result<(), Box<dyn Error>>;
    fn read_parquet(path: &str) -> Result<Self, Box<dyn Error>>;
}

#[cfg(feature="parquet")]
impl SimpleParquet for DataFrame {
    fn write_parquet(&self, path: &str) -> Result<(), Box<dyn Error>> {
        WithParquet::write_parquet(self, path, CompressionOptions::Uncompressed)
    }

    fn read_parquet(path: &str) -> Result<Self, Box<dyn Error>> {
        WithParquet::read_parquet(path)
    }
}
