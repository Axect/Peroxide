use num_complex::Complex;
use crate::traits::math::{Vector, Normed, Norm, InnerProduct};

impl Vector for Complex<f64> {
    type Scalar = Self;

    fn add_vec<'a, 'b>(&'a self, rhs: &'b Self) -> Self {
        self + rhs
    }

    fn sub_vec<'a, 'b>(&'a self, rhs: &'b Self) -> Self {
        self - rhs
    }

    fn mul_scalar(&self, rhs: Self::Scalar) -> Self {
        self * rhs
    }
}

impl Normed for Complex<f64> {
    type UnsignedScalar = f64;
    fn norm(&self, kind: Norm) -> Self::UnsignedScalar {
        match kind {
            Norm::L1 => self.l1_norm(),
            Norm::L2 => Complex::<f64>::norm(*self),
            _ => panic!("No more norms form complex"),
        }
    }

    fn normalize(&self, kind: Norm) -> Self
    where
            Self: Sized {
        let n = self.norm(kind);
        self / n
    }
}

impl InnerProduct for Complex<f64> {
    fn dot(&self, rhs: &Self) -> Self::Scalar {
        self.conj() * rhs
    }
}
