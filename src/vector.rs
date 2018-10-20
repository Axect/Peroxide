use std::convert;
use std::fmt;

#[derive(Clone)]
pub struct Vector {
    data: Vec<f64>,
    len: usize,
}

/// Vector Constructor
pub trait CreateVector<T: convert::Into<f64> + Clone> {
    fn from_vec(v: Vec<T>) -> Vector;
    fn new(start: T, end: T, step: usize) -> Vector;
}

impl<T> CreateVector<T> for Vector where T: convert::Into<f64> + Clone {
    /// from Vec to Vector
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let v = Vector::from_vec(vec![1,2,3,4]);
    /// println!("{}", v); // [1,2,3,4]
    /// ```
    fn from_vec(v: Vec<T>) -> Vector {
        Vector {
            data: v.clone().into_iter().map(|x| x.into()).collect::<Vec<f64>>(),
            len: v.len(),
        }
    }
    /// Vector Constructor from range
    fn new(start: T, end: T, step: usize) -> Vector {
        let s: f64 = start.into();
        let e: f64 = end.into();
        let step: f64 = step as f64;

        assert!(e > s);

        let factor: f64 = (e - s) / step;
        let l: usize = factor as usize + 1;
        let mut v: Vec<f64> = Vec::new();

        for i in 0 .. l {
            v.push(s + step * (i as f64));
        }

        Vector::from_vec(v)
    }
}

impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self.data)
    }
}

impl PartialEq for Vector {
    fn eq(&self, other: &Vector) -> bool {
        self.data == other.data
    }
}

/// R like seq function
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = seq(1,10,1);
/// println!("{}", a); // [1, 2, ..., 10]
/// ```
pub fn seq<T>(start: T, end: T, step: usize) -> Vector where T: convert::Into<f64> + Clone{
    Vector::new(start, end, step)
}
