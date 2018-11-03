extern crate rand;
pub use self::rand::prelude::*;

#[allow(unused_imports)]
use matrix::*;
#[allow(unused_imports)]
use vector::*;

/// MATLAB like zeros
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = zeros!(4);
/// assert_eq!(a, c!(0,0,0,0));
/// ```
#[macro_export]
macro_rules! zeros {
    ( $n:expr ) => {
        vec![0f64; $n]
    };
}

/// MATLAB like rand
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
