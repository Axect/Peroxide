use crate::structure::matrix::{Matrix, Shape};

/// Pointer for col or row
pub trait RawMatrix {
    unsafe fn row_ptr(&self, idx: usize) -> Vec<*const f64>;
    unsafe fn col_ptr(&self, idx: usize) -> Vec<*const f64>;
}

impl RawMatrix for Matrix {
    /// Row pointer
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// let a: Matrix = MATLAB::new("1 2;3 4");
    /// unsafe {
    ///     let b = a.row_ptr(1);
    ///     let b_vec = ptr_to_vec(&b);
    ///     assert_eq!(b_vec, vec![3f64, 4f64]);
    /// }
    /// ```
    unsafe fn row_ptr(&self, idx: usize) -> Vec<*const f64> {
        assert!(idx < self.col, "Index out of range");
        match self.shape {
            Shape::Row => {
                let mut v: Vec<*const f64> = Vec::with_capacity(self.col);
                v.set_len(self.col);
                let start_idx = idx * self.col;
                let p = self.ptr();
                for (i, j) in (start_idx .. start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Col => {
                let mut v: Vec<*const f64> = Vec::with_capacity(self.col);
                v.set_len(self.col);
                let p = self.ptr();
                for (i, elem) in v.iter_mut().enumerate() {
                    *elem = p.add(idx + i * self.row);
                }
                v
            }
        }
    }

    unsafe fn col_ptr(&self, idx: usize) -> Vec<*const f64> {
        assert!(idx < self.col, "Index out of range");
        match self.shape {
            Shape::Col => {
                let mut v: Vec<*const f64> = Vec::with_capacity(self.row);
                v.set_len(self.row);
                let start_idx = idx * self.row;
                let p = self.ptr();
                for (i, j) in (start_idx .. start_idx + v.len()).enumerate() {
                    v[i] = p.add(j);
                }
                v
            }
            Shape::Row => {
                let mut v: Vec<*const f64> = Vec::with_capacity(self.row);
                v.set_len(self.row);
                let p = self.ptr();
                for (i, elem) in v.iter_mut().enumerate() {
                    *elem = p.add(idx + i * self.col);
                }
                v
            }
        }
    }
}
