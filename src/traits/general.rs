/// Some algorithms for Vector
pub trait Algorithm {
    type Scalar;

    fn rank(&self) -> Vec<usize>;
    fn sign(&self) -> Self::Scalar;
    fn arg_max(&self) -> usize;
    fn arg_min(&self) -> usize;
    fn max(&self) -> Self::Scalar;
    fn min(&self) -> Self::Scalar;
    fn swap_with_perm(&mut self, p: &Vec<(usize, usize)>);
}
