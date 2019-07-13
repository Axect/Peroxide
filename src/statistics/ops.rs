/// Factorial
///
/// # Type
/// : `usize -> usize`
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// assert_eq!(factorial(5), 120);
/// ```
pub fn factorial(n: usize) -> usize {
    let mut p = 1usize;
    for i in 1..(n + 1) {
        p *= i;
    }
    p
}

/// Permutation
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// assert_eq!(P(5,3), 60);
/// ```
#[allow(non_snake_case)]
pub fn P(n: usize, r: usize) -> usize {
    let mut p = 1usize;
    for i in 0..r {
        p *= n - i;
    }
    p
}

/// Combination
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// assert_eq!(C(10, 9), 10);
/// ```
#[allow(non_snake_case)]
pub fn C(n: usize, r: usize) -> usize {
    if r > n / 2 {
        return C(n, n - r);
    }

    P(n, r) / factorial(r)
}

/// Combination with Repetition
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// assert_eq!(H(5,3), C(7,3));
/// ```
#[allow(non_snake_case)]
pub fn H(n: usize, r: usize) -> usize {
    C(n + r - 1, r)
}
