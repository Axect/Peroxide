use std::convert;
use std::f64::{MIN};

pub type Vector = Vec<f64>;

pub trait FPVector {
    fn fmap<F>(&self, f: F) -> Vector where F: Fn(f64) -> f64;
    fn reduce<F, T>(&self, init: T, f: F) -> f64
        where F: Fn(f64, f64) -> f64,
              T: convert::Into<f64>;
    fn zip_with<F>(&self, f: F, other: &Vector) -> Vector
        where F: Fn(f64, f64) -> f64;
}

impl FPVector for Vector {
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
    fn fmap<F>(&self, f: F) -> Vector where F: Fn(f64) -> f64 {
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
        where F: Fn(f64, f64) -> f64,
              T: convert::Into<f64> {
        self.clone().into_iter().fold(
            init.into(),
            |x,y| f(x,y),
        )
    }

    fn zip_with<F>(&self, f: F, other: &Vector) -> Vector
        where F: Fn(f64, f64) -> f64 {
        self.into_iter().zip(other).map(|(x,y)| f(*x,*y)).collect::<Vector>()
    }
}

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
        let idx = (1 .. (l+1)).map(|x| x as usize).collect::<Vec<usize>>();

        let mut vec_tup = self.clone().into_iter().zip(idx.clone()).collect::<Vec<(f64, usize)>>();
        vec_tup.sort_by(|x,y| x.0.partial_cmp(&y.0).unwrap().reverse());
        let indices = vec_tup.into_iter().map(|(_,y)| y).collect::<Vec<usize>>();
        idx.into_iter().map(|x| indices.clone().into_iter().position(|t| t == x).unwrap()).collect::<Vec<usize>>()
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

        for i in 0 .. (l-1) {
            for j in 0 .. (l - 1 - i) {
                if v[j] > v[j+1] {
                    sgn *= -1f64;
                    let (a, b) = (v[j], v[j+1]);
                    v[j] = b;
                    v[j+1] = a;
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
        let m = self.clone().into_iter().fold(MIN, |x,y| x.max(y));
        v.into_iter().position(|x| x == m).unwrap()
    }
}

pub trait VecOps {
    type Output;
    fn add(&self, other: &Self) -> Self::Output;
}

impl VecOps for Vector {
    type Output = Vector;
    fn add(&self, other: &Vector) -> Vector {
        self.zip_with(|x, y| x + y, other)
    }
}