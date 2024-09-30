use crate::structure::matrix::{matrix, Matrix, Shape};
use crate::traits::fp::FPVector;
use crate::util::non_macro::zeros_shape;
use std::ops::{Add, Div, Mul, Sub};

/// Syntactic sugar for Vector operations
pub trait VecOps: Sized + FPVector
where
    Self::Scalar: Copy
        + Clone
        + Add<Self::Scalar, Output = Self::Scalar>
        + Sub<Self::Scalar, Output = Self::Scalar>
        + Mul<Self::Scalar, Output = Self::Scalar>
        + Div<Self::Scalar, Output = Self::Scalar>,
{
    //type Scalar;
    fn add_v(&self, v: &Self) -> Self {
        self.zip_with(|x, y| x + y, v)
    }
    fn sub_v(&self, v: &Self) -> Self {
        self.zip_with(|x, y| x - y, v)
    }
    fn mul_v(&self, v: &Self) -> Self {
        self.zip_with(|x, y| x * y, v)
    }
    fn div_v(&self, v: &Self) -> Self {
        self.zip_with(|x, y| x / y, v)
    }
    fn add_s(&self, s: Self::Scalar) -> Self {
        self.fmap(|x| x + s)
    }
    fn sub_s(&self, s: Self::Scalar) -> Self {
        self.fmap(|x| x - s)
    }
    fn mul_s(&self, s: Self::Scalar) -> Self {
        self.fmap(|x| x * s)
    }
    fn div_s(&self, s: Self::Scalar) -> Self {
        self.fmap(|x| x / s)
    }
}

pub trait Scalable {
    type Vec;
    fn reshape(&self, size: (usize, usize), shape: Shape) -> Matrix;
    fn add_col(&self, v: &Self::Vec) -> Matrix;
    fn add_row(&self, v: &Self::Vec) -> Matrix;
}

pub trait ScalableMut {
    type Vec;
    fn reshape_mut(&mut self, size: (usize, usize), shape: Shape);
    fn add_col_mut(&mut self, v: &Self::Vec);
    fn add_row_mut(&mut self, v: &Self::Vec);
}

pub trait ConvToMat {
    fn to_col(&self) -> Matrix;
    fn to_row(&self) -> Matrix;
}

// =============================================================================
// Implementations
// =============================================================================

impl VecOps for Vec<f64> {}
//    /// Vector + Vector
//    ///
//    /// ```
//    /// #[macro_use]
//    /// extern crate peroxide;
//    /// use peroxide::fuga::*;
//    ///
//    /// fn main() {
//    ///     let a = c!(1,2,3,4,5);
//    ///     let b = c!(5,4,3,2,1);
//    ///     assert_eq!(a.add_v(&b), c!(6,6,6,6,6));
//    /// }
//    /// ```
//    fn add_v(&self, v: &Self) -> Self {
//        self.add_vec(&v)
//    }
//
//    /// Vector - Vector
//    ///
//    /// ```
//    /// #[macro_use]
//    /// extern crate peroxide;
//    /// use peroxide::fuga::*;
//    ///
//    /// fn main() {
//    ///     let a = c!(1,2,3,4,5);
//    ///     let b = c!(5,4,3,2,1);
//    ///     assert_eq!(a.sub_v(&b), c!(-4, -2, 0, 2, 4));
//    /// }
//    /// ```
//    fn sub_v(&self, v: &Self) -> Self {
//        self.zip_with(|x, y| x - y, v)
//    }
//
//    /// Vector * Vector
//    ///
//    /// ```
//    /// #[macro_use]
//    /// extern crate peroxide;
//    /// use peroxide::fuga::*;
//    ///
//    /// fn main() {
//    ///     let a = c!(1,2,3,4,5);
//    ///     let b = c!(5,4,3,2,1);
//    ///     assert_eq!(a.mul_v(&b), c!(5, 8, 9, 8, 5));
//    /// }
//    /// ```
//    fn mul_v(&self, v: &Self) -> Self {
//        self.zip_with(|x, y| x * y, v)
//    }
//
//    /// Vector / Vector
//    ///
//    /// ```
//    /// #[macro_use]
//    /// extern crate peroxide;
//    /// use peroxide::fuga::*;
//    ///
//    /// fn main() {
//    ///     let a = c!(2, 4, 6, 8, 10);
//    ///     let b = c!(2, 2, 2, 2, 2);
//    ///     assert_eq!(a.div_v(&b), c!(1,2,3,4,5));
//    /// }
//    /// ```
//    fn div_v(&self, v: &Self) -> Self {
//        self.zip_with(|x, y| x / y, v)
//    }
//
//    /// Vector + Scalar
//    ///
//    /// ```
//    /// #[macro_use]
//    /// extern crate peroxide;
//    /// use peroxide::fuga::*;
//    ///
//    /// fn main() {
//    ///     let a = c!(1,2,3,4,5);
//    ///     let b = 1f64;
//    ///     assert_eq!(a.add_s(b), c!(2,3,4,5,6));
//    /// }
//    /// ```
//    fn add_s(&self, s: Self::Scalar) -> Self {
//        self.fmap(|x| x + s)
//    }
//
//    /// Vector - Scalar
//    ///
//    /// ```
//    /// #[macro_use]
//    /// extern crate peroxide;
//    /// use peroxide::fuga::*;
//    ///
//    /// fn main() {
//    ///     let a = c!(1,2,3,4,5);
//    ///     let b = 1f64;
//    ///     assert_eq!(a.sub_s(b), c!(0,1,2,3,4));
//    /// }
//    /// ```
//    fn sub_s(&self, s: Self::Scalar) -> Self {
//        self.fmap(|x| x - s)
//    }
//
//    /// Vector * Scalar
//    ///
//    /// ```
//    /// #[macro_use]
//    /// extern crate peroxide;
//    /// use peroxide::fuga::*;
//    ///
//    /// fn main() {
//    ///     let a = c!(1,2,3,4,5);
//    ///     let b = 2f64;
//    ///     assert_eq!(a.mul_s(b), c!(2,4,6,8,10));
//    /// }
//    /// ```
//    fn mul_s(&self, s: Self::Scalar) -> Self {
//        self.mul_scalar(s)
//    }
//
//    /// Vector / Scalar
//    ///
//    /// ```
//    /// #[macro_use]
//    /// extern crate peroxide;
//    /// use peroxide::fuga::*;
//    ///
//    /// fn main() {
//    ///     let a = c!(2,4,6,8,10);
//    ///     let b = 2f64;
//    ///     assert_eq!(a.div_s(b), c!(1,2,3,4,5));
//    /// }
//    /// ```
//    fn div_s(&self, s: Self::Scalar) -> Self {
//        self.fmap(|x| x / s)
//    }

impl Scalable for Vec<f64> {
    type Vec = Self;

    /// Vector to Matrix
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3,4,5,6);
    ///     let b1 = a.reshape((3,2), Row);
    ///     let b2 = a.reshape((3,2), Col);
    ///     assert_eq!(b1, ml_matrix("1 2;3 4;5 6"));
    ///     assert_eq!(b2, ml_matrix("1 4;2 5;3 6"));
    /// }
    /// ```
    fn reshape(&self, (r, c): (usize, usize), shape: Shape) -> Matrix {
        assert_eq!(self.len(), r * c);
        let mut m = zeros_shape(r, c, shape);
        m.data = self[..].to_vec();
        m
    }

    /// Vector + Vector = Matrix
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3);
    ///     let b = c!(4,5,6);
    ///     let c1 = a.add_row(&b);
    ///     let c2 = a.add_col(&b);
    ///     assert_eq!(c1, ml_matrix("1 2 3;4 5 6"));
    ///     assert_eq!(c2, ml_matrix("1 4;2 5;3 6"));
    /// }
    /// ```
    fn add_col(&self, v: &Self::Vec) -> Matrix {
        assert_eq!(self.len(), v.len());
        let mut x = self[..].to_vec();
        x.extend_from_slice(&v[..]);
        x.reshape((self.len(), 2), Shape::Col)
    }

    /// Vector + Vector = Matrix
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = c!(1,2,3);
    ///     let b = c!(4,5,6);
    ///     let c1 = a.add_row(&b);
    ///     let c2 = a.add_col(&b);
    ///     assert_eq!(c1, ml_matrix("1 2 3;4 5 6"));
    ///     assert_eq!(c2, ml_matrix("1 4;2 5;3 6"));
    /// }
    /// ```
    fn add_row(&self, v: &Self::Vec) -> Matrix {
        assert_eq!(self.len(), v.len());
        let mut x = self[..].to_vec();
        x.extend_from_slice(&v[..]);
        x.reshape((2, self.len()), Shape::Row)
    }
}

/// Modify Matrix
impl Scalable for Matrix {
    type Vec = Vec<f64>;

    /// Resize matrix
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = ml_matrix("1 2 3;4 5 6"); // ml_matrix has shape `Col`
    ///     let b1 = a.reshape((3, 2), Row);
    ///     let b2 = a.reshape((3, 2), Col);
    ///     assert_eq!(b1, ml_matrix("1 2;3 4;5 6"));
    ///     assert_eq!(b2, ml_matrix("1 4;2 5;3 6"));
    /// }
    /// ```
    fn reshape(&self, (r, c): (usize, usize), shape: Shape) -> Matrix {
        assert_eq!(self.row * self.col, r * c);
        let mut m = zeros_shape(r, c, shape);
        m.data = self.data[..].to_vec();
        m
    }

    /// Add column
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = ml_matrix("1 2 3;4 5 6");
    ///     let b = c!(7,8);
    ///     let c = a.add_col(&b);
    ///     assert_eq!(c, ml_matrix("1 2 3 7;4 5 6 8"));
    /// }
    /// ```
    fn add_col(&self, v: &Self::Vec) -> Matrix {
        assert_eq!(self.row, v.len());
        match self.shape {
            Shape::Col => {
                let mut m = self.clone();
                m.data.extend_from_slice(&v[..]);
                m.col += 1;
                m
            }
            Shape::Row => {
                let mut m = self.change_shape();
                m.data.extend_from_slice(&v[..]);
                m.col += 1;
                m
            }
        }
    }

    /// Add row
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = ml_matrix("1 2 3;4 5 6");
    ///     let b = c!(7,8,9);
    ///     let c = a.add_row(&b);
    ///     assert_eq!(c, ml_matrix("1 2 3;4 5 6;7 8 9"));
    /// }
    /// ```
    fn add_row(&self, v: &Self::Vec) -> Matrix {
        assert_eq!(self.col, v.len());
        match self.shape {
            Shape::Row => {
                let mut m = self.clone();
                m.data.extend_from_slice(&v[..]);
                m.row += 1;
                m
            }
            Shape::Col => {
                let mut m = self.change_shape();
                m.data.extend_from_slice(&v[..]);
                m.row += 1;
                m
            }
        }
    }
}

impl ScalableMut for Matrix {
    type Vec = Vec<f64>;

    /// Resize matrix (Mutable)
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let mut a = ml_matrix("1 2 3;4 5 6"); // ml_matrix has shape `Row`
    ///     a.reshape_mut((3, 2), Row);
    ///     assert_eq!(a, ml_matrix("1 2;3 4;5 6"));
    ///     a.reshape_mut((3, 2), Col);
    ///     assert_eq!(a, ml_matrix("1 4;2 5;3 6"));
    /// }
    /// ```
    fn reshape_mut(&mut self, (r, c): (usize, usize), shape: Shape) {
        assert_eq!(self.row * self.col, r * c);
        self.row = r;
        self.col = c;
        self.shape = shape;
    }

    /// Add column (Mutable)
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let mut a = ml_matrix("1 2 3;4 5 6");
    ///     let b = c!(7,8);
    ///     a.add_col_mut(&b);
    ///     assert_eq!(a, ml_matrix("1 2 3 7;4 5 6 8"));
    /// }
    /// ```
    fn add_col_mut(&mut self, v: &Self::Vec) {
        assert_eq!(self.row, v.len());
        match self.shape {
            Shape::Col => {
                self.data.extend_from_slice(&v[..]);
                self.col += 1;
            }
            Shape::Row => {
                self.change_shape_mut();
                self.data.extend_from_slice(&v[..]);
                self.col += 1;
            }
        }
    }

    /// Add row (Mutable)
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let mut a = ml_matrix("1 2 3;4 5 6");
    ///     let b = c!(7,8,9);
    ///     a.add_row_mut(&b);
    ///     assert_eq!(a, ml_matrix("1 2 3;4 5 6;7 8 9"));
    /// }
    /// ```
    fn add_row_mut(&mut self, v: &Self::Vec) {
        assert_eq!(self.col, v.len());
        match self.shape {
            Shape::Row => {
                self.data.extend_from_slice(&v[..]);
                self.row += 1;
            }
            Shape::Col => {
                self.change_shape_mut();
                self.data.extend_from_slice(&v[..]);
                self.row += 1;
            }
        }
    }
}

impl ConvToMat for Vec<f64> {
    fn to_col(&self) -> Matrix {
        matrix(self.clone(), self.len(), 1, Shape::Col)
    }

    fn to_row(&self) -> Matrix {
        matrix(self.clone(), 1, self.len(), Shape::Row)
    }
}
