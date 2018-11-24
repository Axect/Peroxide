pub trait PowOps {
    type Output;
    fn pow(&self, n: usize) -> Self::Output;
    fn powf(&self, f: f64) -> Self::Output;
}

pub trait TrigOps {
    type Output;
    fn sin(&self) -> Self::Output;
    fn cos(&self) -> Self::Output;
    fn tan(&self) -> Self::Output;
}

pub trait ExpLogOps {
    type Output;
    fn exp(&self) -> Self::Output;
    fn ln(&self) -> Self::Output;
}
