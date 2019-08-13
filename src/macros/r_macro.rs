extern crate rand;
pub use self::rand::prelude::*;

#[allow(unused_imports)]
use statistics::stat::*;
#[allow(unused_imports)]
use structure::matrix::*;
#[allow(unused_imports)]
use structure::vector::*;

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
    ( $start:expr, $end:expr, $step:expr ) => {{
        let s = $start as f64;
        let e = $end as f64;
        let step = $step as f64;

        assert!(e > s);

        let factor: f64 = (e - s) / step;
        let l: usize = factor as usize + 1;
        let mut v: Vec<f64> = Vec::new();

        for i in 0..l {
            v.push(s + step * (i as f64));
        }
        v
    }};
    ( $start:expr; $end:expr; $step:expr ) => {
        seq!($start, $end, $step)
    };
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
///
/// let m = matrix!(0, 2, 2, Row);
/// assert_eq!(m, zeros!(2, 2));
/// ```
#[macro_export]
macro_rules! matrix {
    ( $start:expr;$end:expr;$step:expr, $row:expr, $col:expr, $shape:expr ) => {{
        matrix(seq!($start, $end, $step), $row, $col, $shape)
    }};

    ( $value:expr, $row:expr, $col:expr, $shape:expr ) => {{
        matrix(vec![$value as f64; $row * $col], $row, $col, $shape)
    }};
}

/// R like cbind
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Col);
/// let b = matrix(c!(5,6),2, 1, Col);
/// let c = matrix(c!(7,8),2,1,Col);
/// assert_eq!(cbind!(a.clone(),b.clone()), matrix!(1;6;1,2,3,Col));
/// assert_eq!(cbind!(a,b,c), matrix!(1;8;1, 2, 4, Col));
/// ```
#[macro_export]
macro_rules! cbind {
    // Two
    ( $x:expr, $y:expr ) => {
        {
            let mut temp = $x;
            if temp.shape != Col {
                temp = temp.change_shape();
            }

            let mut v: Vec<f64> = temp.data;
            let mut c: usize = temp.col;
            let r: usize = temp.row;

            assert_eq!(r, $y.row);
            v.extend(&$y.data.clone());
            c += &$y.col;
            matrix(v, r, c, Col)
        }
    };

    // Multi
    ( $x0:expr, $( $x: expr ),* ) => {
        {
            let mut temp0 = $x0;
            if temp0.shape != Col {
                temp0 = temp0.change_shape();
            }
            let mut v: Vec<f64> = temp0.data;
            let mut c: usize = temp0.col;
            let r: usize = temp0.row;
            $(
                let mut temp = $x;
                if temp.shape != Col {
                    temp = temp.change_shape();
                }
                // Must equal row
                assert_eq!(r, temp.row);
                // Add column
                c += temp.col;
                v.extend(&temp.data.clone());
            )*
            matrix(v, r, c, Col)
        }
    };
}

/// R like rbind
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// let b = matrix(c!(5,6),1, 2, Row);
/// let c = matrix(c!(7,8),1, 2, Row);
/// assert_eq!(rbind!(a.clone(),b.clone()), matrix!(1;6;1, 3, 2, Row));
/// assert_eq!(rbind!(a,b,c), matrix!(1;8;1, 4, 2, Row));
/// ```
#[macro_export]
macro_rules! rbind {
    // Two
    ( $x:expr, $y:expr ) => {
        {
            let mut temp = $x;
            if temp.shape != Row {
                temp = temp.change_shape();
            }

            let mut temp2 = $y;
            if temp2.shape != Row {
                temp2 = temp2.change_shape();
            }

            let mut v: Vec<f64> = temp.data;
            let c: usize = temp.col;
            let mut r: usize = temp.row;

            assert_eq!(c, temp2.col);
            v.extend(&temp2.data.clone());
            r += temp2.row;
            matrix(v, r, c, Row)
        }
    };

    // Multi
    ( $x0:expr, $( $x: expr ),* ) => {
        {
            let mut temp0 = $x0;
            if temp0.shape != Row {
                temp0 = temp0.change_shape();
            }
            let mut v: Vec<f64> = temp0.data;
            let c: usize = temp0.col;
            let mut r: usize = temp0.row;
            $(
                let mut temp = $x;
                if temp.shape != Row {
                    temp = temp.change_shape();
                }
                // Must equal row
                assert_eq!(c, temp.col);
                // Add column
                r += temp.row;
                v.extend(&temp.data.clone());
            )*
            matrix(v, r, c, Row)
        }
    };
}

/// R like random uniform
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = runif!(5, -1, 1);
/// println!("{:?}", a);
///
/// let b = runif!(5); // same as runif!(5,0,1)
/// println!("{:?}", b);
/// ```
#[macro_export]
macro_rules! runif {
    ( $x0:expr, $start:expr, $end:expr ) => {{
        let n: usize = $x0;
        let mut v = vec![0f64; n];

        let mut rng = thread_rng();

        for i in 0..n {
            v[i] = rng.gen_range($start as f64, $end as f64);
        }
        v
    }};

    ( $x0:expr ) => {{
        runif!($x0, 0, 1)
    }};
}

/// R like lm
///
/// # Caution
/// Should use single token to use this macro
/// * `lm(y ~ x)` can be accepted
/// * `lm(c!(1,2,3,4) ~ c!(1,2,3,4))` is not allowed!
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
/// use peroxide::*;
///
/// let a = matrix!(1;5;1, 5, 1, Col);
/// let b = matrix(c!(3.7, 4.2, 4.9, 5.7, 6.0), 5, 1, Col);
///
/// assert_eq!(lm!(b ~ a), matrix(c!(3.0700, 0.6100), 2, 1, Col));
/// ```
#[macro_export]
macro_rules! lm {
    ( $y:tt ~ $x:tt ) => {
        lm(&$x, &$y)
    };
}
