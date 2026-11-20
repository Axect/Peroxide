/// Factorial
///
/// # Type
/// : `usize -> usize`
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// assert_eq!(factorial(5), 120);
/// ```
///
/// # Panics
/// Panics if the result overflows `usize` (`n > 20` on 64-bit).
pub fn factorial(n: usize) -> usize {
    let mut p = 1usize;
    for i in 1..(n + 1) {
        p = p.checked_mul(i).expect("factorial: overflow");
    }
    p
}

/// Double Factorial
///
/// # Type
/// : `usize -> usize`
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// assert_eq!(double_factorial(7), 105);
/// ```
pub fn double_factorial(n: usize) -> usize {
    let mut s = 1usize;
    let mut n = n;
    while n >= 2 {
        s = s.checked_mul(n).expect("double_factorial: overflow");
        n -= 2;
    }
    s
}

/// Permutation
///
/// Returns 0 when `r > n` (the falling factorial contains a zero factor).
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// assert_eq!(P(5,3), 60);
/// assert_eq!(P(3,5), 0);
/// ```
///
/// # Panics
/// Panics if the result overflows `usize`.
#[allow(non_snake_case)]
pub fn P(n: usize, r: usize) -> usize {
    if r > n {
        return 0;
    }
    let mut p = 1usize;
    for i in 0..r {
        p = p.checked_mul(n - i).expect("P: overflow");
    }
    p
}

/// Combination (binomial coefficient)
///
/// Returns 0 when `r > n`.
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// assert_eq!(C(10, 9), 10);
/// assert_eq!(C(30, 15), 155117520);
/// assert_eq!(C(5, 7), 0);
/// ```
///
/// # Panics
/// Panics if an intermediate product overflows `usize`. The multiplicative
/// formula keeps intermediates no larger than `r * C(n, r)`, so every
/// `C(n, r)` with `n <= 62` is exact on 64-bit.
#[allow(non_snake_case)]
pub fn C(n: usize, r: usize) -> usize {
    if r > n {
        return 0;
    }
    let r = r.min(n - r);
    let mut c = 1usize;
    for i in 1..=r {
        // c * (n - r + i) is divisible by i: c = C(n - r + i - 1, i - 1) here
        c = c.checked_mul(n - r + i).expect("C: overflow") / i;
    }
    c
}

/// Combination with Repetition
///
/// # Usage
///
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// assert_eq!(H(5,3), C(7,3));
/// ```
#[allow(non_snake_case)]
pub fn H(n: usize, r: usize) -> usize {
    C(n + r - 1, r)
}
