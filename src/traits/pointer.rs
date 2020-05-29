use traits::math::Vector;
use std::ops::{Add, Deref, Sub, Mul, Div};
use traits::fp::FPVector;
use Dual;

// =============================================================================
// Redox Structure
// =============================================================================
#[derive(Debug)]
pub struct Redox<T: Vector> {
    data: Box<T>
}

impl<T: Vector> Deref for Redox<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl Redox<Vec<f64>> {
    pub fn from_vec(vec: Vec<f64>) -> Self {
        Self {
            data: Box::new(vec)
        }
    }

    pub fn red(self) -> Vec<f64> {
        (*self).to_vec()
    }
}

impl Redox<Vec<Dual>> {
    pub fn from_vec(vec: Vec<Dual>) -> Self {
        Self {
            data: Box::new(vec)
        }
    }

    pub fn red(self) -> Vec<Dual> {
        (*self).to_vec()
    }
}

// =============================================================================
// Oxide trait
// =============================================================================
pub trait Oxide: Vector {
    fn ox(self) -> Redox<Self> where Self: Sized;
}

// =============================================================================
// Reference Arithmetics
// =============================================================================

impl<T: Vector> Add<Redox<T>> for Redox<T> {
    type Output = Self;

    fn add(self, rhs: Redox<T>) -> Self::Output {
        Redox {
            data: Box::new(self.add_vec(&rhs.data))
        }
    }
}

impl<T: Vector> Sub<Redox<T>> for Redox<T> {
    type Output = Self;

    fn sub(self, rhs: Redox<T>) -> Self::Output {
        Redox {
            data: Box::new(self.add_vec(&rhs.data))
        }
    }
}

impl<T: Vector + FPVector> Mul<Redox<T>> for Redox<T> 
where
    <T as FPVector>::Scalar: Mul<Output=<T as FPVector>::Scalar>
{
    type Output = Self;

    fn mul(self, rhs: Redox<T>) -> Self::Output {
        Redox {
            data: Box::new(self.zip_with(|x, y| x * y, &rhs.data))
        }
    }
}

impl<T: Vector + FPVector> Div<Redox<T>> for Redox<T> 
where
    <T as FPVector>::Scalar: Div<Output=<T as FPVector>::Scalar>
{
    type Output = Self;

    fn div(self, rhs: Redox<T>) -> Self::Output {
        Redox {
            data: Box::new(self.zip_with(|x, y| x / y, &rhs.data))
        }
    }
}

impl<T: Vector + FPVector> Add<f64> for Redox<T>
where
    <T as FPVector>::Scalar: Add<f64, Output=<T as FPVector>::Scalar> 
{
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        Redox {
            data: Box::new(self.fmap(|x| x + rhs))
        }
    }
}

impl<T: Vector + FPVector> Sub<f64> for Redox<T>
where
    <T as FPVector>::Scalar: Sub<f64, Output=<T as FPVector>::Scalar> 
{
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        Redox {
            data: Box::new(self.fmap(|x| x - rhs))
        }
    }
}

impl<T: Vector + FPVector> Mul<f64> for Redox<T>
where
    <T as FPVector>::Scalar: Mul<f64, Output=<T as FPVector>::Scalar> 
{
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Redox {
            data: Box::new(self.fmap(|x| x * rhs))
        }
    }
}

impl<T: Vector + FPVector> Div<f64> for Redox<T>
where
    <T as FPVector>::Scalar: Div<f64, Output=<T as FPVector>::Scalar> 
{
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Redox {
            data: Box::new(self.fmap(|x| x / rhs))
        }
    }
}