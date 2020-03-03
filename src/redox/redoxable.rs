use std::ops::{Add, Deref, Index, IndexMut, Sub, Neg, Mul, Div};
use structure::vector::VecOps;
use std::net::Shutdown::Read;
use PowOps;

/// Smart pointer of Vector
pub struct RedoxVector {
    pub data: Vec<f64>
}

impl RedoxVector {
    pub fn new(n: usize) -> Self {
        RedoxVector {
            data: vec![f64::default(); n]
        }
    }

    pub fn from_vec(v: Vec<f64>) -> Self {
        RedoxVector {
            data: v
        }
    }
}

impl Deref for RedoxVector {
    type Target = Vec<f64>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl Index<usize> for RedoxVector {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        self.data.index(index)
    }
}

impl IndexMut<usize> for RedoxVector {
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        self.data.index_mut(index)
    }
}

impl Neg for RedoxVector {
    type Output = RedoxVector;

    fn neg(self) -> Self::Output {
        RedoxVector::from_vec(self.data.into_iter().map(|x| -x).collect())
    }
}

impl Add<f64> for RedoxVector {
    type Output = RedoxVector;

    fn add(self, rhs: f64) -> Self::Output {
        RedoxVector::from_vec(self.data.s_add(rhs))
    }
}

impl Sub<f64> for RedoxVector {
    type Output = RedoxVector;

    fn sub(self, rhs: f64) -> Self::Output {
        RedoxVector::from_vec(self.data.s_sub(rhs))
    }
}

impl Mul<f64> for RedoxVector {
    type Output = RedoxVector;

    fn mul(self, rhs: f64) -> Self::Output {
        RedoxVector::from_vec(self.data.s_mul(rhs))
    }
}

impl Div<f64> for RedoxVector {
    type Output = RedoxVector;

    fn div(self, rhs: f64) -> Self::Output {
        RedoxVector::from_vec(self.data.s_div(rhs))
    }
}

impl Add<RedoxVector> for f64 {
    type Output = RedoxVector;

    fn add(self, rhs: RedoxVector) -> Self::Output {
        rhs.add(self)
    }
}

impl Sub<RedoxVector> for f64 {
    type Output = RedoxVector;

    fn sub(self, rhs: RedoxVector) -> Self::Output {
        RedoxVector::from_vec(rhs.data.into_iter().map(|x| self - x).collect())
    }
}

impl Mul<RedoxVector> for f64 {
    type Output = RedoxVector;

    fn mul(self, rhs: RedoxVector) -> Self::Output {
        rhs.mul(self)
    }
}

impl Add<RedoxVector> for RedoxVector {
    type Output = RedoxVector;

    fn add(self, rhs: Self) -> Self::Output {
        RedoxVector::from_vec(self.data.add(&rhs.data))
    }
}

impl<'a, 'b> Add<&'b RedoxVector> for &'a RedoxVector {
    type Output = RedoxVector;

    fn add(self, rhs: &'b RedoxVector) -> Self::Output {
        RedoxVector::from_vec(self.data.add(&rhs.data))
    }
}

impl Sub<RedoxVector> for RedoxVector {
    type Output = RedoxVector;

    fn sub(self, rhs: Self) -> Self::Output {
        RedoxVector::from_vec(self.data.sub(&rhs.data))
    }
}

impl<'a, 'b> Sub<&'b RedoxVector> for &'a RedoxVector {
    type Output = RedoxVector;

    fn sub(self, rhs: &'b RedoxVector) -> Self::Output {
        RedoxVector::from_vec(self.data.sub(&rhs.data))
    }
}

/// Dot product
impl Mul<RedoxVector> for RedoxVector {
    type Output = f64;

    fn mul(self, rhs: Self) -> Self::Output {
        self.data.dot(&rhs.data)
    }
}

impl<'a, 'b> Mul<&'b RedoxVector> for &'a RedoxVector {
    type Output = f64;

    fn mul(self, rhs: &'b RedoxVector) -> Self::Output {
        self.data.dot(&rhs.data)
    }
}

/// Redox wrap for Vector
pub trait Redoxable {
    fn redox(self) -> RedoxVector;
}

impl Redoxable for Vec<f64> {
    fn redox(self) -> RedoxVector {
        RedoxVector::from_vec(self)
    }
}