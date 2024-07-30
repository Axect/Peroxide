//! Useful missing tools

use std::ops::Range;

// =============================================================================
// Fundamental Utils
// =============================================================================

/// Near equal
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// assert!(nearly_eq(1.0/3.0 * 3.0, 1));
/// ```
pub fn nearly_eq<S, T>(x: S, y: T) -> bool
where
    S: Into<f64>,
    T: Into<f64>,
{
    let mut b: bool = false;
    let e = 1e-7;
    let p: f64 = x.into().abs();
    let q: f64 = y.into().abs();
    if (p - q).abs() < e || (p - q).abs() / (p + q).min(f64::MAX) < e {
        b = true;
    }
    b
}

#[allow(unused_comparisons)]
pub fn tab(s: &str, space: usize) -> String {
    let l = s.len();
    let mut m: String = String::new();
    let fs = format!("{}{}", " ".repeat(space - l), s);
    m.push_str(&fs);
    m
}

pub fn quot_rem(x: usize, y: usize) -> (i32, i32) {
    ((x / y) as i32, (x % y) as i32)
}

// =============================================================================
// Choose Utils
// =============================================================================
pub fn choose_shorter_string(x1: String, x2: String) -> String {
    if x1.len() > x2.len() {
        x2
    } else {
        x1
    }
}

pub fn choose_shorter_vec(x1: &Vec<f64>, x2: &Vec<f64>) -> Vec<f64> {
    if x1.len() > x2.len() {
        x2.clone()
    } else {
        x1.clone()
    }
}

pub fn choose_longer_vec(x1: &Vec<f64>, x2: &Vec<f64>) -> Vec<f64> {
    if x1.len() <= x2.len() {
        x2.clone()
    } else {
        x1.clone()
    }
}

pub fn max<T>(v: Vec<T>) -> T
where
    T: PartialOrd + Copy + Clone,
{
    let l = v.len();
    if l == 1 {
        v[0]
    } else {
        let mut t = if v[0] >= v[1] { v[0] } else { v[1] };
        for i in 2..v.len() {
            if v[i] > t {
                t = v[i];
            }
        }
        t
    }
}

pub fn min<T>(v: Vec<T>) -> T
where
    T: PartialOrd + Copy + Clone,
{
    let l = v.len();
    if l == 1 {
        v[0]
    } else {
        let mut t = if v[0] <= v[1] { v[0] } else { v[1] };
        for i in 2..v.len() {
            if v[i] < t {
                t = v[i];
            }
        }
        t
    }
}

/// Signum function
pub fn sgn(x: usize) -> f64 {
    if x % 2 == 0 {
        1f64
    } else {
        -1f64
    }
}

/// Vector compare
pub fn eq_vec(x: &[f64], y: &[f64], tol: f64) -> bool {
    x.iter().zip(y.iter()).all(|(x, y)| (x - y).abs() <= tol)
}

// =============================================================================
// Vec of Tuples
// =============================================================================
/// Auto-zip
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let a = vec![1, 2, 3];
/// let a_zipped = auto_zip(&a);
/// assert_eq!(a_zipped, vec![(1, 2), (2, 3)]);
/// ```
pub fn auto_zip<T: Clone>(x: &Vec<T>) -> Vec<(T, T)> {
    let x_head = x[0..x.len() - 1].to_vec();
    let x_tail = x[1..x.len()].to_vec();
    x_head.into_iter().zip(x_tail).collect()
}

/// Find the index of interval of x
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// let x = vec![
///     (0, 5),
///     (5, 7),
///     (7, 10),
///     (10, 15),
///     (15, 20),
///     (20, 30)
/// ];
///
/// assert_eq!(find_interval(&x, 11), 3);
/// ```
pub fn find_interval<T: PartialOrd + PartialEq>(sorted_intervals: &Vec<(T, T)>, x: T) -> usize {
    let mut i = 0;
    let mut j = sorted_intervals.len() - 1;

    // Check range
    assert!(
        x >= sorted_intervals[0].0,
        "x is smaller than the smallest interval"
    );
    assert!(
        x <= sorted_intervals[sorted_intervals.len() - 1].1,
        "x is larger than the largest interval"
    );

    while i <= j {
        let mid = (i + j) / 2;
        if x < sorted_intervals[mid].0 {
            j = mid - 1;
        } else if x > sorted_intervals[mid].1 {
            i = mid + 1;
        } else {
            return mid;
        }
    }
    i
}

/// Generate Range of Intervals
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use std::ops::Range;
///
/// let x = vec![1, 2, 3, 4];
/// let r = gen_range(&x);
///
/// let answer = vec![
///     Range { start: 1, end: 2 },
///     Range { start: 2, end: 3 },
///     Range { start: 3, end: 4 },
/// ];
///
/// assert_eq!(r, answer);
/// ```
pub fn gen_range<T: Clone + PartialOrd>(x: &[T]) -> Vec<Range<T>> {
    let mut r = Vec::new();
    for i in 0..x.len() - 1 {
        r.push(Range {
            start: x[i].clone(),
            end: x[i + 1].clone(),
        });
    }
    r
}

/// Generate and Zip Range of Intervals
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::fuga::*;
/// use std::ops::Range;
///
/// let x: Vec<i32> = vec![1, 2, 3, 4];
/// let y: Vec<i32> = x.iter().map(|&t| t.pow(2)).collect();
/// let r = zip_range(&x, &y);
///
/// let answer = vec![
///     (Range { start: 1, end: 2 }, 1),
///     (Range { start: 2, end: 3 }, 4),
///     (Range { start: 3, end: 4 }, 9),
/// ];
///
/// assert_eq!(r, answer);
/// ```
pub fn zip_range<T: Clone + PartialOrd, U: Clone>(x: &[T], y: &[U]) -> Vec<(Range<T>, U)> {
    y[0..x.len() - 1]
        .iter()
        .enumerate()
        .map(|(i, yi)| {
            (
                Range {
                    start: x[i].clone(),
                    end: x[i + 1].clone(),
                },
                yi.clone(),
            )
        })
        .collect()
}
