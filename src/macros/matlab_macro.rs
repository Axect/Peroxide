/// MATLAB like zeros - zero matrix
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = zeros!(4);
///     assert_eq!(a, c!(0,0,0,0));
///
///     let b = zeros!(3, 2);
///     assert_eq!(b, matrix(c!(0,0,0,0,0,0), 3, 2, Row));
/// }
/// ```
#[macro_export]
macro_rules! zeros {
    ( $n:expr ) => {
        vec![0f64; $n]
    };

    ( $r:expr, $c:expr ) => {{
        let (r, c) = ($r, $c);
        matrix(vec![0f64; r * c], r, c, Row)
    }};
}

/// MATLAB like rand - random matrix
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = rand!(2, 2);
///     println!("{}", a);  // 2 x 2 random matrix (0 ~ 1)
/// }
/// ```
#[macro_export]
macro_rules! rand {
    () => {{
        let mut rng = rand::rng();
        rng.random_range(0f64..=1f64)
    }};

    ( $m:expr, $n:expr ) => {{
        let r = $m;
        let c = $n;
        let mut rng = rand::rng();
        let mut m = matrix(vec![0f64; r * c], r, c, Row);
        for i in 0..r {
            for j in 0..c {
                m[(i, j)] = rng.random_range(0f64..=1f64);
            }
        }
        m
    }};
}

/// MATLAB like eye - identity matrix
///
/// # Examples
///
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let i = eye!(2);
///     assert_eq!(i, matrix(c!(1,0,0,1), 2, 2, Row));
/// }
/// ```
#[macro_export]
macro_rules! eye {
    ( $n:expr ) => {{
        let n = $n;
        let mut m = matrix(vec![0f64; n * n], n, n, Row);
        for i in 0..n {
            m[(i, i)] = 1f64;
        }
        m
    }};
}

/// MATLAB like linspace
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = linspace!(1, 10, 10);
///     assert_eq!(a, seq!(1,10,1));
/// }
/// ```
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = linspace!(10, 1000);
///     assert_eq!(a, seq!(10,1000,10));
/// }
/// ```
#[macro_export]
macro_rules! linspace {
    ( $start:expr, $end:expr, $length: expr) => {{
        let step = ($end - $start) as f64 / ($length as f64 - 1f64);
        seq!($start, $end, step)
    }};

    ( $start:expr, $end:expr ) => {{
        let step = ($end - $start) as f64 / (99f64);
        seq!($start, $end, step)
    }};
}
