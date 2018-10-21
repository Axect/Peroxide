/// R like concatenate
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

/// R like seq macro (Not precise)
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// assert_eq!(seq!(1,10,1), c!(1,2,3,4,5,6,7,8,9,10));
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
}
