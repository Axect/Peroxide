extern crate rand;
pub use self::rand::prelude::*;

#[allow(unused_imports)]
use matrix::*;
#[allow(unused_imports)]
use vector::*;
use r_macro::*;

/// MATLAB like zeros - zero matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = zeros!(4);
/// assert_eq!(a, c!(0,0,0,0));
///
/// let b = zeros!(3, 2);
/// assert_eq!(b, matrix(c!(0,0,0,0,0,0), 3, 2, Row));
/// ```
#[macro_export]
macro_rules! zeros {
    ( $n:expr ) => {
        vec![0f64; $n]
    };

    ( $r:expr, $c:expr ) => {
        {
            let (r, c) = ($r, $c);
            matrix(vec![0f64; r*c], r, c, Row)
        }
    };
}

/// MATLAB like rand - random matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
/// 
/// let a = rand!(2, 2);
/// println!("{}", a); // 2 x 2 random matrix (0 ~ 1)
/// ```
#[macro_export]
macro_rules! rand {
    ( $m:expr, $n:expr ) => {
        {   
            let r = $m;
            let c = $n;
            let mut rng = thread_rng();
            let mut m = matrix(vec![0f64; r * c], r, c, Row);
            for i in 0 .. r {
                for j in 0 .. c {
                    m[(i, j)] = rng.gen_range(0f64, 1f64);
                }
            }
            m
        }
    };
}

/// MATLAB like eye - identity matrix
///
/// # Examples
///
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let i = eye!(2);
/// assert_eq!(i, matrix(c!(1,0,0,1), 2, 2, Row));
/// ```
#[macro_export]
macro_rules! eye {
    ( $n:expr ) => {
        {
            let n = $n;
            let mut m = matrix(vec![0f64; n * n], n, n, Row);
            for i in 0 .. n {
                m[(i, i)] = 1f64;
            }
            m
        }
    };
}

/// MATLAB like linspace
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = linspace(1, 10, 10);
/// assert_eq!(a, seq!(1,10,1));
/// ```
macro_rules! linspace {
    ( $start:expr, $end:expr, $length: expr) => {
        {
            let step = ($end - $start) / ($length - 1f64);
            seq!($start, $end, step)
        }
    };
}