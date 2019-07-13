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
