use operation::extra_ops::PowOps;
use std::convert;
use std::f64::MIN;

pub type Vector = Vec<f64>;

/// Functional Programming tools for Vector
pub trait FPVector {
    type Scalar;

    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> Self::Scalar;
    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
        T: convert::Into<Self::Scalar>;
    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar;
    fn filter<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> bool;
    fn take(&self, n: usize) -> Self;
    fn skip(&self, n: usize) -> Self;
}

impl FPVector for Vector {
    type Scalar = f64;

    /// fmap for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// assert_eq!(a.fmap(|x| x*2f64), seq!(2,10,2));
    /// ```
    fn fmap<F>(&self, f: F) -> Vector
    where
        F: Fn(f64) -> f64,
    {
        self.clone().into_iter().map(|x| f(x)).collect::<Vector>()
    }

    /// reduce for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = seq!(1,100,1);
    /// assert_eq!(a.reduce(0, |x,y| x + y), 5050f64);
    /// ```
    fn reduce<F, T>(&self, init: T, f: F) -> f64
    where
        F: Fn(f64, f64) -> f64,
        T: convert::Into<f64>,
    {
        self.clone().into_iter().fold(init.into(), |x, y| f(x, y))
    }

    fn zip_with<F>(&self, f: F, other: &Vector) -> Vector
    where
        F: Fn(f64, f64) -> f64,
    {
        self.into_iter()
            .zip(other)
            .map(|(x, y)| f(*x, *y))
            .collect::<Vector>()
    }

    /// Filter for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// let b = a.filter(|x| x > 3.);
    /// assert_eq!(b, c!(4,5));
    /// ```
    fn filter<F>(&self, f: F) -> Vector
    where
        F: Fn(f64) -> bool,
    {
        self.clone()
            .into_iter()
            .filter(|x| f(*x))
            .collect::<Vector>()
    }

    /// Take for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// let b = a.take(3);
    /// assert_eq!(b, c!(1,2,3));
    /// ```
    fn take(&self, n: usize) -> Vector {
        let mut v = vec![0f64; n];
        for i in 0..n {
            v[i] = self[i];
        }
        return v;
    }

    /// Skip for Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,2,3,4,5);
    /// let b = a.skip(3);
    /// assert_eq!(b, c!(4,5));
    /// ```
    fn skip(&self, n: usize) -> Vector {
        let l = self.len();
        let mut v = vec![0f64; l - n];
        for i in n..l {
            v[i - n] = self[i];
        }
        return v;
    }
}

/// Some algorithms for Vector
pub trait Algorithm {
    fn rank(&self) -> Vec<usize>;
    fn sign(&self) -> f64;
    fn arg_max(&self) -> usize;
}

impl Algorithm for Vector {
    /// Assign rank
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let v = c!(7, 5, 9, 2, 8);
    /// assert_eq!(v.rank(), vec![2,3,0,4,1]);
    /// ```
    fn rank(&self) -> Vec<usize> {
        let l = self.len();
        let idx = (1..(l + 1)).map(|x| x as usize).collect::<Vec<usize>>();

        let mut vec_tup = self
            .clone()
            .into_iter()
            .zip(idx.clone())
            .collect::<Vec<(f64, usize)>>();
        vec_tup.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap().reverse());
        let indices = vec_tup.into_iter().map(|(_, y)| y).collect::<Vec<usize>>();
        idx.into_iter()
            .map(|x| indices.clone().into_iter().position(|t| t == x).unwrap())
            .collect::<Vec<usize>>()
    }

    /// Sign of Permutation
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = c!(1,0,2);
    /// let b = c!(1,2,0);
    /// let c = c!(0,1,2);
    ///
    /// assert_eq!((a.sign(), b.sign(), c.sign()), (-1f64, 1f64, 1f64));
    /// ```
    fn sign(&self) -> f64 {
        let l = self.len();
        let mut v = self.clone();
        let mut sgn = 1f64;

        for i in 0..(l - 1) {
            for j in 0..(l - 1 - i) {
                if v[j] > v[j + 1] {
                    sgn *= -1f64;
                    let (a, b) = (v[j], v[j + 1]);
                    v[j] = b;
                    v[j + 1] = a;
                }
            }
        }
        return sgn;
    }

    /// arg max
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let v = c!(1,3,2,4,3,7);
    /// assert_eq!(v.arg_max(),5);
    ///
    /// let v2 = c!(1,3,2,5,6,6);
    /// assert_eq!(v2.arg_max(),4);
    /// ```
    fn arg_max(&self) -> usize {
        let v = self.clone();
        let m = self.clone().into_iter().fold(MIN, |x, y| x.max(y));
        v.into_iter().position(|x| x == m).unwrap()
    }
}

/// Convenient Vector Operation trait
pub trait VecOps {
    type Scalar;
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn div(&self, other: &Self) -> Self;
    fn dot(&self, other: &Self) -> Self::Scalar;
    fn norm(&self) -> Self::Scalar;
}

/// Convenient Vector Operations (No Clone, No Copy)
impl VecOps for Vector {
    type Scalar = f64;

    /// Addition
    fn add(&self, other: &Vector) -> Vector {
        self.zip_with(|x, y| x + y, other)
    }

    /// Subtraction
    fn sub(&self, other: &Vector) -> Vector {
        self.zip_with(|x, y| x - y, other)
    }

    /// Multiplication
    fn mul(&self, other: &Vector) -> Vector {
        self.zip_with(|x, y| x * y, other)
    }

    /// Division
    fn div(&self, other: &Vector) -> Vector {
        self.zip_with(|x, y| x / y, other)
    }

    /// Dot product
    fn dot(&self, other: &Vector) -> f64 {
        self.mul(other).reduce(0, |x, y| x + y)
    }

    /// Norm
    fn norm(&self) -> f64 {
        self.dot(&self).sqrt()
    }
}

/// Power operation for Vector
impl PowOps for Vector {
    /// Power usize
    fn powi(&self, n: i32) -> Self {
        self.powf(n as f64)
    }

    /// Power float
    fn powf(&self, f: f64) -> Self {
        self.fmap(|x| x.powf(f))
    }

    /// Sqrt
    fn sqrt(&self) -> Self {
        self.powf(0.5)
    }
}
