/// Some algorithms for Vector
pub trait Algorithm {
    fn rank(&self) -> Vec<usize>;
    fn sign(&self) -> f64;
    fn arg_max(&self) -> usize;
    fn swap_with_perm(&mut self, p: &Vec<(usize, usize)>);
}
