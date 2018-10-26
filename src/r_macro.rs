#[allow(unused_imports)]
use matrix::*;
#[allow(unused_imports)]
use vector::*;

/// R like concatenate (Type: Vec\<f64\>)
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = c![1,2,3,4];
/// let b = c![5,6,7,8];
/// let c = c![a; b];
/// println!("{:?}", a); // [1,2,3,4]
/// println!("{:?}", b); // [5,6,7,8]
/// println!("{:?}", c); // [1,2,3,4,5,6,7,8]
/// ```
#[macro_export]
macro_rules! c {
    ( $( $x:expr ),* ) => {
        {
            let mut v: Vec<f64> = Vec::new();
            let mut l: usize = 0;
            $(
                v.push($x as f64);
                l += 1;
            )*
            v
        }
    };
    ( $( $x:expr );* ) => {
        {
            let mut v: Vec<f64> = Vec::new();
            $(
                v.extend(&$x);
            )*
            v
        }
    }
}

/// R like seq macro
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// assert_eq!(seq!(1,10,1), c!(1,2,3,4,5,6,7,8,9,10));
/// assert_eq!(seq!(1,10,1), seq!(1;10;1));
/// ```
#[macro_export]
macro_rules! seq {
    ( $start:expr, $end:expr, $step:expr ) => {
        {
            let s = $start as f64;
            let e = $end as f64;
            let step = $step as f64;

            assert!(e > s);

            let factor: f64 = (e - s) / step;
            let l: usize = factor as usize + 1;
            let mut v: Vec<f64> = Vec::new();

            for i in 0 .. l {
                v.push(s + step * (i as f64));
            }
            v
        }
    };
    ( $start:expr; $end:expr; $step:expr ) => {
        seq!($start, $end, $step)
    }
}

/// More R like Matrix constructor (Macro)
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row); // start;end;step
/// let b = matrix(c!(1,2,3,4), 2, 2, Row);
/// let c = matrix(vec![1,2,3,4], 2, 2, Row); // Normal function
/// assert!(a == b && b == c);
/// ```
#[macro_export]
macro_rules! matrix {
    ( $start:expr;$end:expr;$step:expr,$row:expr,$col:expr,$shape:expr ) => {
        {
            matrix(
                seq!($start,$end,$step),
                $row,
                $col,
                $shape
            )
        }
    };
}

// Python like
/// zeros - like numpy
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
