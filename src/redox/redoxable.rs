use std::ops::{Add, Deref, Index, IndexMut};

pub struct RedoxVector<T> {
    pub data: Vec<T>
}

impl<T: Default + Clone> RedoxVector<T> {
    pub fn new(n: usize) -> Self {
        RedoxVector {
            data: vec![T::default(); n]
        }
    }

    pub fn from_vec(v: Vec<T>) -> Self {
        RedoxVector {
            data: v
        }
    }
}

impl<T> Deref for RedoxVector<T> {
    type Target = Vec<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T> Index<usize> for RedoxVector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        self.data.index(index)
    }
}

impl<T> IndexMut<usize> for RedoxVector<T> {
    fn index_mut(&mut self, index: usize) -> &mut T {
        self.data.index_mut(index)
    }
}

impl<T> Add<RedoxVector<T>> for RedoxVector<T>
where T: Default + Add<Output=T> + Copy + Clone
{
    type Output = RedoxVector<T>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut rv: RedoxVector<T> = RedoxVector::new((*self).len());
        for i in 0 .. rv.len() {
            rv[i] = self[i] + rhs[i];
        }
        rv
    }
}

impl<'a, 'b, T> Add<&'b RedoxVector<T>> for &'a RedoxVector<T>
where T: Default + Add<Output=T> + Copy + Clone
{
    type Output = RedoxVector<T>;

    fn add(self, rhs: &'b RedoxVector<T>) -> Self::Output {
        let mut rv: RedoxVector<T> = RedoxVector::new((*self).len());
        for i in 0 .. rv.len() {
            rv[i] = self[i] + rhs[i];
        }
        rv
    }
}