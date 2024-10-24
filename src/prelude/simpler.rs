use crate::numerical::{
    eigen,
    eigen::{Eigen, EigenMethod::Jacobi},
    integral,
    integral::Integral::G7K15R,
    spline,
    spline::{CubicHermiteSpline, SlopeMethod::Quadratic},
};
#[cfg(feature = "parquet")]
use crate::structure::dataframe::{DataFrame, WithParquet};
use crate::structure::matrix::Matrix;
use crate::structure::polynomial;
use crate::traits::math::{Norm, Normed};
use crate::traits::matrix::{MatrixTrait, LinearAlgebra, PQLU, WAZD, QR, Form, SolveKind};
#[cfg(feature = "parquet")]
use arrow2::io::parquet::write::CompressionOptions;
#[cfg(feature = "parquet")]
use std::error::Error;

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
pub trait SimplerLinearAlgebra<M:MatrixTrait> {
    fn back_subs(&self, b: &[f64]) -> Vec<f64>;
    fn forward_subs(&self, b: &[f64]) -> Vec<f64>;
    fn lu(&self) -> PQLU<M>;
    fn waz_diag(&self) -> Option<WAZD<M>>;
    fn waz(&self) -> Option<WAZD<M>>;
    fn qr(&self) -> QR<M>;
    #[cfg(feature = "O3")]
    fn cholesky(&self) -> M;
    fn rref(&self) -> M;
    fn det(&self) -> f64;
    fn block(&self) -> (M, M, M, M);
    fn inv(&self) -> M;
    fn pseudo_inv(&self) -> M;
    fn solve(&self, b: &[f64]) -> Vec<f64>;
    fn solve_mat(&self, m: &M) -> M;
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

impl SimplerLinearAlgebra<Matrix> for Matrix {
    fn back_subs(&self, b: &[f64]) -> Vec<f64> {
        LinearAlgebra::back_subs(self, b)
    }

    fn forward_subs(&self, b: &[f64]) -> Vec<f64> {
        LinearAlgebra::forward_subs(self, b)
    }

    fn lu(&self) -> PQLU<Matrix> {
        LinearAlgebra::lu(self)
    }

    fn waz_diag(&self) -> Option<WAZD<Matrix>> {
        LinearAlgebra::waz(self, Form::Diagonal)
    }

    fn waz(&self) -> Option<WAZD<Matrix>> {
        LinearAlgebra::waz(self, Form::Identity)
    }

    fn qr(&self) -> QR<Matrix> {
        LinearAlgebra::qr(self)
    }

    #[cfg(feature = "O3")]
    fn cholesky(&self) -> Matrix {
        LinearAlgebra::cholesky(self, UPLO::Lower)
    }

    fn rref(&self) -> Matrix {
        LinearAlgebra::rref(self)
    }

    fn det(&self) -> f64 {
        LinearAlgebra::det(self)
    }

    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix) {
        LinearAlgebra::block(self)
    }

    fn inv(&self) -> Matrix {
        LinearAlgebra::inv(self)
    }

    fn pseudo_inv(&self) -> Matrix {
        LinearAlgebra::pseudo_inv(self)
    }

    fn solve(&self, b: &[f64]) -> Vec<f64> {
        LinearAlgebra::solve(self, b, SolveKind::LU)
    }

    fn solve_mat(&self, m: &Matrix) -> Matrix {
        LinearAlgebra::solve_mat(self, m, SolveKind::LU)
    }

    fn is_symmetric(&self) -> bool {
        LinearAlgebra::is_symmetric(self)
    }
}

/// Simple solve
#[allow(non_snake_case)]
pub fn solve(A: &Matrix, m: &Matrix) -> Matrix {
    crate::traits::matrix::solve(A, m, SolveKind::LU)
}

/// Simple Chebyshev Polynomial (First Kind)
pub fn chebyshev_polynomial(n: usize) -> polynomial::Polynomial {
    polynomial::chebyshev_polynomial(n, polynomial::SpecialKind::First)
}

pub fn cubic_hermite_spline(node_x: &[f64], node_y: &[f64]) -> anyhow::Result<CubicHermiteSpline> {
    spline::cubic_hermite_spline(node_x, node_y, Quadratic)
}

use crate::special::function::{
    lambert_w0 as lambert_w0_flex, lambert_wm1 as lambert_wm1_flex, LambertWAccuracyMode,
};

/// The principal branch of the Lambert W function, W_0(`z`).
///
/// Returns [`NAN`](f64::NAN) if the given input is smaller than -1/e (≈ -0.36787944117144233).
///
/// Accurate to 50 bits.
///
/// Wrapper of `lambert_w_0` function of `lambert_w` crate.
///
/// # Reference
///
/// [Toshio Fukushima, Precise and fast computation of Lambert W function by piecewise minimax rational function approximation with variable transformation](https://www.researchgate.net/publication/346309410_Precise_and_fast_computation_of_Lambert_W_function_by_piecewise_minimax_rational_function_approximation_with_variable_transformation)
pub fn lambert_w0(z: f64) -> f64 {
    lambert_w0_flex(z, LambertWAccuracyMode::Precise)
}

/// The secondary branch of the Lambert W function, W_-1(`z`).
///
/// Returns [`NAN`](f64::NAN) if the given input is positive or smaller than -1/e (≈ -0.36787944117144233).
///
/// Accurate to 50 bits.
///
/// Wrapper of `lambert_w_m1` function of `lambert_w` crate.
///
/// # Reference
///
/// [Toshio Fukushima, Precise and fast computation of Lambert W function by piecewise minimax rational function approximation with variable transformation](https://www.researchgate.net/publication/346309410_Precise_and_fast_computation_of_Lambert_W_function_by_piecewise_minimax_rational_function_approximation_with_variable_transformation)
pub fn lambert_wm1(z: f64) -> f64 {
    lambert_wm1_flex(z, LambertWAccuracyMode::Precise)
}

/// Simple handle parquet
#[cfg(feature = "parquet")]
pub trait SimpleParquet: Sized {
    fn write_parquet(&self, path: &str) -> Result<(), Box<dyn Error>>;
    fn read_parquet(path: &str) -> Result<Self, Box<dyn Error>>;
}

#[cfg(feature = "parquet")]
impl SimpleParquet for DataFrame {
    fn write_parquet(&self, path: &str) -> Result<(), Box<dyn Error>> {
        WithParquet::write_parquet(self, path, CompressionOptions::Uncompressed)
    }

    fn read_parquet(path: &str) -> Result<Self, Box<dyn Error>> {
        WithParquet::read_parquet(path)
    }
}
