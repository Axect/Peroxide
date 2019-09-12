use std::ops::{Add, Deref, Index, IndexMut};
use structure::vector::VecOps;

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