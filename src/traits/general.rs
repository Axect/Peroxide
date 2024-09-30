/// Some algorithms for Vector
pub trait Algorithm {
    fn rank(&self) -> Vec<usize>;
    fn sign(&self) -> f64;
    fn arg_max(&self) -> usize;
    fn arg_min(&self) -> usize;
    fn max(&self) -> f64;
    fn min(&self) -> f64;
    fn swap_with_perm(&mut self, p: &Vec<(usize, usize)>);
}

/// Some algorithms for Vector in Parallel
pub trait ParallelAlgorithm {
    fn par_rank(&self) -> Vec<usize>;
    fn par_arg_max(&self) -> usize;
    fn par_arg_min(&self) -> usize;
    fn par_max(&self) -> f64;
    fn par_min(&self) -> f64;

    // Not implemented in parallel:
    // fn par_sign(&self) -> f64;
    // fn par_swap_with_perm(&mut self, p: &Vec<(usize, usize)>);
}
