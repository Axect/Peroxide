use num_complex::Complex;
use crate::traits::fp::FPVector;
use crate::traits::math::{Vector, Normed, Norm, InnerProduct};
use crate::traits::sugar::VecOps;

impl Vector for Complex<f64> {
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

impl Normed for Complex<f64> {
    type UnsignedScalar = f64;
    fn norm(&self, kind: Norm) -> Self::UnsignedScalar {
        match kind {
            Norm::L1 => self.l1_norm(),
            Norm::L2 => Complex::<f64>::norm(*self),
            _ => unimplemented!(),
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

impl FPVector for Vec<Complex<f64>> {
    type Scalar = Complex<f64>;

    fn fmap<F>(&self, f: F) -> Self
    where
            F: Fn(Self::Scalar) -> Self::Scalar {
        self.iter().map(|&x| f(x)).collect()
    }

    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
            F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar {
        self.iter().zip(other.iter()).map(|(&x, &y)| f(x,y)).collect()
    }

    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
            F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
            T: Into<Self::Scalar> {
        self.iter().fold(init.into(), |x, &y| f(x,y))
    }

    fn filter<F>(&self, f: F) -> Self
    where
            F: Fn(Self::Scalar) -> bool {
        self.iter().filter(|&x| f(*x)).cloned().collect()
    }

    fn take(&self, n: usize) -> Self {
        self.iter().take(n).cloned().collect()
    }

    fn skip(&self, n: usize) -> Self {
        self.iter().skip(n).cloned().collect()
    }

    fn sum(&self) -> Self::Scalar {
        self.iter().sum()
    }

    fn prod(&self) -> Self::Scalar {
        self.iter().product()
    }
}

impl Vector for Vec<Complex<f64>> {
    type Scalar = Complex<f64>;

    fn add_vec(&self, rhs: &Self) -> Self {
        self.zip_with(|x, y| x + y, rhs)
    }

    fn sub_vec(&self, rhs: &Self) -> Self {
        self.zip_with(|x, y| x - y, rhs)
    }

    fn mul_scalar(&self, rhs: Self::Scalar) -> Self {
        self.fmap(|x| x * rhs)
    }
}

impl Normed for Vec<Complex<f64>> {
    type UnsignedScalar = f64;

    fn norm(&self, kind: Norm) -> Self::UnsignedScalar {
        match kind {
            Norm::L1 => self.iter().map(|x| Complex::<f64>::norm(*x).abs()).sum(),
            _ => unimplemented!()
        }
    }

    fn normalize(&self, kind: Norm) -> Self
    where
            Self: Sized {
        unimplemented!()
    }
}

impl InnerProduct for Vec<Complex<f64>> {
    fn dot(&self, rhs: &Self) -> Self::Scalar {
        self.zip_with(|x, y| x.conj() * y, rhs).sum()
    }
}

impl VecOps for Vec<Complex<f64>> {}
