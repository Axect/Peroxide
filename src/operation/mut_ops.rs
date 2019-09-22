use ::{Matrix, Shape};
use util::low_level::swap_vec_ptr;

pub trait MutFP {
    type Scalar;
    fn mut_map<F>(&mut self, f: F)
    where
        F: Fn(Self::Scalar) -> Self::Scalar;
    fn mut_zip_with<F>(&mut self, f: F, other: &Self)
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar;
}

impl MutFP for Vec<f64> {
    type Scalar = f64;

    fn mut_map<F>(&mut self, f: F)
    where
        F: Fn(Self::Scalar) -> Self::Scalar,
    {
        for i in 0..self.len() {
            self[i] = f(self[i]);
        }
    }

    fn mut_zip_with<F>(&mut self, f: F, other: &Self)
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
    {
        for i in 0..self.len() {
            self[i] = f(self[i], other[i]);
        }
    }
}

pub trait MutMatrix {
    unsafe fn col_mut(&mut self, idx: usize) -> Vec<*mut f64>;
    unsafe fn row_mut(&mut self, idx: usize) -> Vec<*mut f64>;
    unsafe fn swap_mut(&mut self, idx1: usize, idx2: usize, shape: Shape);
}

impl MutMatrix for Matrix {
    unsafe fn col_mut(&mut self, idx: usize) -> Vec<*mut f64> {
        assert!(idx < self.col, "Index out of range");
        match self.shape {
            Shape::Col => {
                let mut v: Vec<*mut f64> = Vec::with_capacity(self.row);
                v.set_len(self.row);
                let start_idx = idx * self.row;
                let p = self.mut_ptr();
                for (i, j) in (start_idx .. start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Row => {
                let mut v: Vec<*mut f64> = Vec::with_capacity(self.row);
                v.set_len(self.row);
                let p = self.mut_ptr();
                for i in 0 .. v.len() {
                    v[i] = p.add(idx + i * self.col);
                }
                v
            }
        }
    }

    unsafe fn row_mut(&mut self, idx: usize) -> Vec<*mut f64> {
        assert!(idx < self.row, "Index out of range");
        match self.shape {
            Shape::Row => {
                let mut v: Vec<*mut f64> = Vec::with_capacity(self.col);
                v.set_len(self.col);
                let start_idx = idx * self.col;
                let p = self.mut_ptr();
                for (i, j) in (start_idx .. start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Col => {
                let mut v: Vec<*mut f64> = Vec::with_capacity(self.col);
                v.set_len(self.col);
                let p = self.mut_ptr();
                for i in 0 .. v.len() {
                    v[i] = p.add(idx + i * self.row);
                }
                v
            }
        }
    }

    unsafe fn swap_mut(&mut self, idx1: usize, idx2: usize, shape: Shape) {
        match shape {
            Shape::Col => {
                swap_vec_ptr(&mut self.col_mut(idx1), &mut self.col_mut(idx2))
            }
            Shape::Row => {
                swap_vec_ptr(&mut self.row_mut(idx1), &mut self.row_mut(idx2))
            }
        }
    }
}
