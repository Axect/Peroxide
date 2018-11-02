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