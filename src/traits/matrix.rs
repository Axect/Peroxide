pub trait MatrixTrait: Sized {
    type Scalar;
    fn ptr(&self) -> *const Self::Scalar; 
    fn mut_ptr(&mut self) -> *mut Self::Scalar;
    fn as_slice(&self) -> &[Self::Scalar];
    fn as_mut_slice(&mut self) -> &mut [Self::Scalar];
    fn change_shape(&self) -> Self;
    fn change_shape_mut(&mut self);
    fn spread(&self) -> String;
    fn col(&self, index: usize) -> Vec<Self::Scalar>;
    fn row(&self, index: usize) -> Vec<Self::Scalar>;
    fn diag(&self) -> Vec<Self::Scalar>;
    fn transpose(&self) -> Self;
    fn t(&self) -> Self { self.transpose() }
    fn subs_col(&mut self, idx: usize, v: &[Self::Scalar]);
    fn subs_row(&mut self, idx: usize, v: &[Self::Scalar]);
    fn from_index<F, G>(f: F, size: (usize, usize)) -> Self
    where
        F: Fn(usize, usize) -> G + Copy,
        G: Into<Self::Scalar>;
    fn to_vec(&self) -> Vec<Vec<Self::Scalar>>;
    fn to_diag(&self) -> Self;
    fn submat(&self, start: (usize, usize), end: (usize, usize)) -> Self;
    fn subs_mat(&mut self, start: (usize, usize), end: (usize, usize), m: &Self);
}
