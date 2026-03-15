//! Taylor mode forward automatic differentiation with const-generic `Jet<N>` type
//!
//! ## Overview
//!
//! This module provides a const-generic `Jet<N>` struct for Taylor-mode forward AD
//! of arbitrary order N. The struct stores **normalized** Taylor coefficients:
//!
//! ```text
//! Jet { value: c_0, deriv: [c_1, c_2, ..., c_N] }
//! ```
//!
//! where $c_k = f^{(k)}(a) / k!$ is the $k$-th normalized Taylor coefficient evaluated
//! at the expansion point $a$. This normalization eliminates binomial coefficients
//! from all arithmetic recurrences.
//!
//! ## Type Aliases
//!
//! * `Dual = Jet<1>` — first-order forward AD (value + first derivative)
//! * `HyperDual = Jet<2>` — second-order forward AD (value + first + second derivative)
//!
//! ## Constructors
//!
//! * `Jet::var(x)` — independent variable at point x (deriv\[0\] = 1)
//! * `Jet::constant(x)` — constant (all derivatives zero)
//! * `Jet::new(value, deriv)` — raw constructor
//! * `ad0(x)` — `Jet<0>` constant (backward compat)
//! * `ad1(x, dx)` — `Jet<1>` with first derivative (backward compat)
//! * `ad2(x, dx, ddx)` — `Jet<2>` with first and second derivatives (backward compat)
//!
//! ## Accessors
//!
//! * `.value()` / `.x()` — $f(a)$
//! * `.dx()` — $f'(a)$
//! * `.ddx()` — $f''(a)$
//! * `.derivative(k)` — $f^{(k)}(a)$ (raw factorial-scaled derivative)
//! * `.taylor_coeff(k)` — normalized Taylor coefficient $c_k$
//!
//! ## Implemented Operations
//!
//! * `Add, Sub, Mul, Div` (Jet op Jet, Jet op f64, f64 op Jet)
//! * `Neg`
//! * `ExpLogOps`: `exp`, `ln`, `log`, `log2`, `log10`
//! * `PowOps`: `powi`, `powf`, `pow`, `sqrt`
//! * `TrigOps`: `sin_cos`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`,
//!              `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`
//!
//! ## Usage
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // First derivative of f(x) = x^2 at x = 2
//!     let x = Jet::<1>::var(2.0);
//!     let y = x.powi(2);
//!     assert_eq!(y.value(), 4.0);
//!     assert_eq!(y.dx(), 4.0);  // f'(2) = 2*2 = 4
//!
//!     // Second derivative using HyperDual
//!     let x2 = HyperDual::new(2.0, [1.0, 0.0]);
//!     let y2 = x2.powi(2);
//!     assert_eq!(y2.value(), 4.0);
//!     assert_eq!(y2.dx(), 4.0);   // f'(2) = 4
//!     assert_eq!(y2.ddx(), 2.0);  // f''(2) = 2
//! }
//! ```
//!
//! ### Higher-order derivatives
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // 5th derivative of x^5 at x = 1
//!     let x = Jet::<5>::var(1.0);
//!     let y = x.powi(5);
//!     assert_eq!(y.derivative(5), 120.0);  // 5! = 120
//! }
//! ```
//!
//! ### Using the `#[ad_function]` macro
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! #[ad_function]
//! fn f(x: f64) -> f64 {
//!     x.sin() + x.powi(2)
//! }
//!
//! fn main() {
//!     // f_grad and f_hess are generated automatically
//!     let grad = f_grad(1.0);   // f'(1) = cos(1) + 2
//!     let hess = f_hess(1.0);   // f''(1) = -sin(1) + 2
//!
//!     assert!((grad - (1.0_f64.cos() + 2.0)).abs() < 1e-10);
//!     assert!((hess - (-1.0_f64.sin() + 2.0)).abs() < 1e-10);
//! }
//! ```
//!
//! ### Generic functions with `Real` trait
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn quadratic<T: Real>(x: T) -> T {
//!     x.powi(2) + x * 3.0 + T::from_f64(1.0)
//! }
//!
//! fn main() {
//!     // Works with both f64 and AD (= Jet<2>)
//!     let val = quadratic(2.0_f64);                // 11.0
//!     let jet = quadratic(AD1(2.0, 1.0));
//!     assert_eq!(val, 11.0);                       // f(2) = 4 + 6 + 1
//!     assert_eq!(jet.value(), 11.0);               // f(2) = 11
//!     assert_eq!(jet.dx(), 7.0);                   // f'(2) = 2*2 + 3
//! }
//! ```
//!
//! ### Jacobian computation
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // Jacobian of f(x,y) = [x - y, x + 2*y] at (1, 1)
//!     let x = vec![1.0, 1.0];
//!     let j = jacobian(f, &x);
//!     j.print();
//!     //       c[0] c[1]
//!     // r[0]     1   -1
//!     // r[1]     1    2
//! }
//!
//! fn f(xs: &Vec<AD>) -> Vec<AD> {
//!     let x = xs[0];
//!     let y = xs[1];
//!     vec![x - y, x + 2.0 * y]
//! }
//! ```
//!
//! ### Backward-compatible constructors
//!
//! ```
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // These work just like the old AD1/AD2 constructors
//!     let a = AD1(2.0, 1.0);   // value=2, f'=1
//!     let b = AD2(4.0, 4.0, 2.0);  // x^2 at x=2
//!
//!     assert_eq!(a.x(), 2.0);
//!     assert_eq!(b.dx(), 4.0);
//!     assert_eq!(b.ddx(), 2.0);
//!
//!     // New constructors (equivalent)
//!     let c = Jet::<1>::var(2.0);  // Same as Dual var at x=2
//!     assert_eq!(c.dx(), 1.0);     // dx/dx = 1 for independent variable
//! }
//! ```
//!
//! ## Accuracy: Jet\<N\> vs Finite Differences
//!
//! `Jet<N>` computes derivatives to **machine precision** because it propagates
//! exact Taylor coefficients through the computation graph. In contrast, finite
//! difference methods suffer from both truncation and cancellation errors that
//! worsen rapidly at higher derivative orders.
//!
//! The plot below compares the relative error of `Jet<N>` against central finite
//! differences ($h = 10^{-4}$) for $f(x) = \sin(x)$ at $x = 1.0$, across derivative orders 1–8:
//!
//! ![Derivative Accuracy](https://raw.githubusercontent.com/Axect/Peroxide/master/example_data/derivative_accuracy.png)
//!
//! `Jet<N>` (blue) stays at $\sim 10^{-15}$ (machine epsilon) for all orders,
//! while finite differences (green) degrade from $\sim 10^{-9}$ at order 1 to $> 10^{0}$ at order 4.
//!
//! ## Taylor Series Convergence
//!
//! Since `Jet<N>` stores normalized Taylor coefficients $c_k = f^{(k)}(a)/k!$,
//! you can directly reconstruct the Taylor polynomial of any function:
//!
//! $$T_N(x) = c_0 + c_1 (x-a) + c_2 (x-a)^2 + \cdots + c_N (x-a)^N$$
//!
//! The plot below shows the Taylor polynomial of $\sin(x)$ around $x = 0$ for
//! increasing truncation orders $N = 1, 3, 5, 7, 9$:
//!
//! ![Taylor Convergence](https://raw.githubusercontent.com/Axect/Peroxide/master/example_data/taylor_convergence.png)
//!
//! As $N$ increases, the Taylor polynomial converges to the exact $\sin(x)$ curve
//! over a wider interval.

use crate::traits::{fp::FPVector, math::Vector, stable::StableFn, sugar::VecOps};
use peroxide_num::{ExpLogOps, PowOps, TrigOps};
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

// =============================================================================
// Jet struct
// =============================================================================

/// Const-generic Taylor-mode forward AD type.
///
/// Stores the value and $N$ normalized Taylor coefficients:
/// - `value` = $f(a) = c_0$
/// - `deriv[k]` = $f^{(k+1)}(a) / (k+1)! = c_{k+1}$
///
/// So `Jet<1>` stores $(c_0, c_1) = (f(a),\, f'(a))$,
/// and `Jet<2>` stores $(c_0, c_1, c_2) = (f(a),\, f'(a),\, f''(a)/2)$.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Jet<const N: usize> {
    value: f64,
    deriv: [f64; N],
}

impl<const N: usize> Jet<N> {
    /// Create a `Jet` from raw value and normalized Taylor coefficient array.
    pub fn new(value: f64, deriv: [f64; N]) -> Self {
        Self { value, deriv }
    }

    /// Create an independent variable jet at point `x`.
    /// Sets `deriv[0] = 1.0` (the 1st normalized coefficient), rest zero.
    ///
    /// # Examples
    /// ```
    /// use peroxide::fuga::*;
    ///
    /// let x = Jet::<2>::var(3.0);
    /// assert_eq!(x.value(), 3.0);
    /// assert_eq!(x.dx(), 1.0);    // dx/dx = 1
    /// assert_eq!(x.ddx(), 0.0);   // d²x/dx² = 0
    /// ```
    pub fn var(x: f64) -> Self {
        let mut deriv = [0.0f64; N];
        if N >= 1 {
            deriv[0] = 1.0;
        }
        Self { value: x, deriv }
    }

    /// Create a constant jet (all derivatives zero).
    pub fn constant(x: f64) -> Self {
        Self {
            value: x,
            deriv: [0.0f64; N],
        }
    }

    /// The function value $f(a)$.
    #[inline]
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Alias for `value()` — backward compatibility.
    #[inline]
    pub fn x(&self) -> f64 {
        self.value
    }

    /// First derivative $f'(a)$.
    ///
    /// # Examples
    /// ```
    /// use peroxide::fuga::*;
    ///
    /// let x = Jet::<1>::var(2.0);
    /// let y = x.powi(3);    // x^3
    /// assert_eq!(y.dx(), 12.0);  // 3*x^2 = 3*4 = 12
    /// ```
    #[inline]
    pub fn dx(&self) -> f64 {
        if N >= 1 {
            self.deriv[0]
        } else {
            0.0
        }
    }

    /// Second derivative $f''(a)$.
    ///
    /// # Examples
    /// ```
    /// use peroxide::fuga::*;
    ///
    /// let x = Jet::<2>::var(2.0);
    /// let y = x.powi(3);     // x^3
    /// assert_eq!(y.ddx(), 12.0);  // 6*x = 6*2 = 12
    /// ```
    #[inline]
    pub fn ddx(&self) -> f64 {
        if N >= 2 {
            self.deriv[1] * 2.0
        } else {
            0.0
        }
    }

    /// Returns $f^{(\mathrm{order})}(a)$, the raw (factorial-scaled) derivative of given order.
    /// - order = 0: $f(a)$
    /// - order = 1: $f'(a)$ = `deriv[0]`
    /// - order = k: $f^{(k)}(a)$ = `deriv[k-1]` $\times\, k!$
    ///
    /// Internally computes `taylor_coeff(k)` $\times\, k!$.
    ///
    /// # Examples
    /// ```
    /// use peroxide::fuga::*;
    ///
    /// let x = Jet::<3>::var(0.0);
    /// let y = x.exp();
    /// // All derivatives of exp at 0 are 1
    /// assert!((y.derivative(0) - 1.0).abs() < 1e-15);
    /// assert!((y.derivative(1) - 1.0).abs() < 1e-15);
    /// assert!((y.derivative(2) - 1.0).abs() < 1e-15);
    /// assert!((y.derivative(3) - 1.0).abs() < 1e-15);
    /// ```
    pub fn derivative(&self, order: usize) -> f64 {
        if order == 0 {
            self.value
        } else if order <= N {
            self.deriv[order - 1] * factorial(order) as f64
        } else {
            0.0
        }
    }

    /// Returns the $k$-th normalized Taylor coefficient $c_k = f^{(k)}(a) / k!$.
    /// - $k = 0$: $c_0 = f(a)$
    /// - $k \ge 1$: $c_k$ = `deriv[k-1]`
    pub fn taylor_coeff(&self, k: usize) -> f64 {
        self.coeff(k)
    }

    /// Internal: get the $k$-th Taylor coefficient $c_k$.
    #[inline]
    fn coeff(&self, k: usize) -> f64 {
        if k == 0 {
            self.value
        } else if k <= N {
            self.deriv[k - 1]
        } else {
            0.0
        }
    }

    /// Internal: set the $k$-th Taylor coefficient $c_k$.
    #[inline]
    fn set_coeff(&mut self, k: usize, v: f64) {
        if k == 0 {
            self.value = v;
        } else if k <= N {
            self.deriv[k - 1] = v;
        }
    }

    /// Internal: create a zero jet.
    #[inline]
    fn zero() -> Self {
        Self {
            value: 0.0,
            deriv: [0.0f64; N],
        }
    }
}

// =============================================================================
// Type aliases
// =============================================================================

/// First-order forward AD: stores value and first derivative.
pub type Dual = Jet<1>;

/// Second-order forward AD: stores value, first derivative, and second derivative $/\, 2!$.
pub type HyperDual = Jet<2>;

// =============================================================================
// Compatibility constructors
// =============================================================================

/// Create a `Jet<0>` constant (zero-order, value only).
#[inline]
pub fn ad0(x: f64) -> Jet<0> {
    Jet { value: x, deriv: [] }
}

/// Create a `Jet<1>` with value and first derivative.
/// `dx` is the raw first derivative $f'(a)$; stored as `deriv[0]` $= dx / 1! = dx$.
///
/// # Arguments
/// * `x` - function value $f(a)$
/// * `dx` - first derivative $f'(a)$
///
/// # Examples
/// ```
/// use peroxide::fuga::*;
///
/// let j = ad1(2.0, 1.0);  // variable x at x=2
/// assert_eq!(j.value(), 2.0);
/// assert_eq!(j.dx(), 1.0);
/// ```
#[inline]
pub fn ad1(x: f64, dx: f64) -> Jet<1> {
    Jet {
        value: x,
        deriv: [dx],
    }
}

/// Create a `Jet<2>` with value, first derivative, and second derivative.
/// `ddx` is the raw second derivative $f''(a)$; stored internally as `deriv[1]` $= f''(a) / 2!$.
///
/// # Arguments
/// * `x` - function value $f(a)$
/// * `dx` - first derivative $f'(a)$
/// * `ddx` - second derivative $f''(a)$
///
/// # Examples
/// ```
/// use peroxide::fuga::*;
///
/// // Represent x^2 at x=2: f=4, f'=4, f''=2
/// let j = ad2(4.0, 4.0, 2.0);
/// assert_eq!(j.value(), 4.0);
/// assert_eq!(j.dx(), 4.0);
/// assert_eq!(j.ddx(), 2.0);
/// ```
#[inline]
pub fn ad2(x: f64, dx: f64, ddx: f64) -> Jet<2> {
    Jet {
        value: x,
        deriv: [dx, ddx / 2.0],
    }
}

// =============================================================================
// Helper: factorial
// =============================================================================

#[inline]
fn factorial(n: usize) -> u64 {
    let mut result = 1u64;
    for i in 2..=(n as u64) {
        result *= i;
    }
    result
}

// =============================================================================
// Display
// =============================================================================

impl<const N: usize> std::fmt::Display for Jet<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Jet({}", self.value)?;
        if N > 0 {
            write!(f, "; ")?;
            for (i, d) in self.deriv.iter().enumerate() {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", d)?;
            }
        }
        write!(f, ")")
    }
}

// =============================================================================
// PartialOrd (compare by value only)
// =============================================================================

impl<const N: usize> PartialOrd for Jet<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.value.partial_cmp(&other.value)
    }
}

// =============================================================================
// From conversions
// =============================================================================

impl<const N: usize> From<f64> for Jet<N> {
    fn from(v: f64) -> Self {
        Self::constant(v)
    }
}

impl<const N: usize> From<Jet<N>> for f64 {
    fn from(j: Jet<N>) -> f64 {
        j.value
    }
}

// =============================================================================
// Index / IndexMut (backward compat: index 0 = value, index k >= 1 = deriv[k-1])
// =============================================================================

impl<const N: usize> Index<usize> for Jet<N> {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        if index == 0 {
            &self.value
        } else if index <= N {
            &self.deriv[index - 1]
        } else {
            panic!("Jet<{}> index {} out of bounds (max index = {})", N, index, N)
        }
    }
}

impl<const N: usize> IndexMut<usize> for Jet<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index == 0 {
            &mut self.value
        } else if index <= N {
            &mut self.deriv[index - 1]
        } else {
            panic!("Jet<{}> index {} out of bounds (max index = {})", N, index, N)
        }
    }
}

// =============================================================================
// Neg
// =============================================================================

impl<const N: usize> Neg for Jet<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut z = self;
        z.value = -z.value;
        for d in z.deriv.iter_mut() {
            *d = -*d;
        }
        z
    }
}

// =============================================================================
// Add, Sub, Mul, Div for Jet<N> op Jet<N>
// =============================================================================

impl<const N: usize> Add<Jet<N>> for Jet<N> {
    type Output = Self;

    fn add(self, rhs: Jet<N>) -> Self::Output {
        let mut z = self;
        z.value += rhs.value;
        for i in 0..N {
            z.deriv[i] += rhs.deriv[i];
        }
        z
    }
}

impl<const N: usize> Sub<Jet<N>> for Jet<N> {
    type Output = Self;

    fn sub(self, rhs: Jet<N>) -> Self::Output {
        let mut z = self;
        z.value -= rhs.value;
        for i in 0..N {
            z.deriv[i] -= rhs.deriv[i];
        }
        z
    }
}

impl<const N: usize> Mul<Jet<N>> for Jet<N> {
    type Output = Self;

    /// Multiplication using normalized Taylor coefficient convolution:
    /// $z_n = \sum_{k=0}^{n} c_k \cdot d_{n-k}$.
    /// No binomial coefficients needed due to normalization convention.
    fn mul(self, rhs: Jet<N>) -> Self::Output {
        let mut z = Self::zero();
        for n in 0..=N {
            let mut s = 0.0f64;
            for k in 0..=n {
                s += self.coeff(k) * rhs.coeff(n - k);
            }
            z.set_coeff(n, s);
        }
        z
    }
}

impl<const N: usize> Div<Jet<N>> for Jet<N> {
    type Output = Self;

    /// Division using normalized Taylor coefficient recurrence:
    /// $z_0 = a_0 / b_0$,
    /// $z_n = \frac{1}{b_0}\left(a_n - \sum_{k=1}^{n} b_k \, z_{n-k}\right)$
    fn div(self, rhs: Jet<N>) -> Self::Output {
        let b0 = rhs.coeff(0);
        let inv_b0 = 1.0 / b0;
        let mut z = Self::zero();
        z.set_coeff(0, self.coeff(0) * inv_b0);
        for n in 1..=N {
            let mut s = 0.0f64;
            for k in 1..=n {
                s += rhs.coeff(k) * z.coeff(n - k);
            }
            z.set_coeff(n, inv_b0 * (self.coeff(n) - s));
        }
        z
    }
}

// =============================================================================
// Scalar arithmetic: Jet<N> op f64
// =============================================================================

impl<const N: usize> Add<f64> for Jet<N> {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        let mut z = self;
        z.value += rhs;
        z
    }
}

impl<const N: usize> Sub<f64> for Jet<N> {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        let mut z = self;
        z.value -= rhs;
        z
    }
}

impl<const N: usize> Mul<f64> for Jet<N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        let mut z = self;
        z.value *= rhs;
        for d in z.deriv.iter_mut() {
            *d *= rhs;
        }
        z
    }
}

impl<const N: usize> Div<f64> for Jet<N> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        let inv = 1.0 / rhs;
        let mut z = self;
        z.value *= inv;
        for d in z.deriv.iter_mut() {
            *d *= inv;
        }
        z
    }
}

// =============================================================================
// Scalar arithmetic: f64 op Jet<N>
// =============================================================================

impl<const N: usize> Add<Jet<N>> for f64 {
    type Output = Jet<N>;

    fn add(self, rhs: Jet<N>) -> Self::Output {
        let mut z = rhs;
        z.value += self;
        z
    }
}

impl<const N: usize> Sub<Jet<N>> for f64 {
    type Output = Jet<N>;

    fn sub(self, rhs: Jet<N>) -> Self::Output {
        let mut z = -rhs;
        z.value += self;
        z
    }
}

impl<const N: usize> Mul<Jet<N>> for f64 {
    type Output = Jet<N>;

    fn mul(self, rhs: Jet<N>) -> Self::Output {
        rhs * self
    }
}

impl<const N: usize> Div<Jet<N>> for f64 {
    type Output = Jet<N>;

    fn div(self, rhs: Jet<N>) -> Self::Output {
        Jet::<N>::constant(self) / rhs
    }
}

// =============================================================================
// ExpLogOps
// =============================================================================

impl<const N: usize> ExpLogOps for Jet<N> {
    type Float = f64;

    /// $\exp(a)$ using the normalized recurrence:
    /// $z_0 = e^{a_0}$,
    /// $z_n = \frac{1}{n}\sum_{k=1}^{n} k\, a_k\, z_{n-k}$
    fn exp(&self) -> Self {
        let mut z = Self::zero();
        z.set_coeff(0, self.coeff(0).exp());
        for n in 1..=N {
            let mut s = 0.0f64;
            for k in 1..=n {
                s += (k as f64) * self.coeff(k) * z.coeff(n - k);
            }
            z.set_coeff(n, s / (n as f64));
        }
        z
    }

    /// $\ln(a)$ using the normalized recurrence:
    /// $z_0 = \ln(a_0)$,
    /// $z_n = \frac{1}{a_0}\left(a_n - \frac{1}{n}\sum_{k=1}^{n-1} k\, z_k\, a_{n-k}\right)$
    fn ln(&self) -> Self {
        let a0 = self.coeff(0);
        let inv_a0 = 1.0 / a0;
        let mut z = Self::zero();
        z.set_coeff(0, a0.ln());
        for n in 1..=N {
            let mut s = 0.0f64;
            for k in 1..n {
                s += (k as f64) * z.coeff(k) * self.coeff(n - k);
            }
            z.set_coeff(n, inv_a0 * (self.coeff(n) - s / (n as f64)));
        }
        z
    }

    fn log(&self, base: f64) -> Self {
        let ln_base = base.ln();
        let z = self.ln();
        let mut result = Self::zero();
        result.set_coeff(0, z.coeff(0) / ln_base);
        for k in 1..=N {
            result.set_coeff(k, z.coeff(k) / ln_base);
        }
        result
    }

    fn log2(&self) -> Self {
        self.log(2.0)
    }

    fn log10(&self) -> Self {
        self.log(10.0)
    }
}

// =============================================================================
// PowOps
// =============================================================================

impl<const N: usize> PowOps for Jet<N> {
    type Float = f64;

    /// Integer power via repeated multiplication.
    fn powi(&self, n: i32) -> Self {
        if n == 0 {
            return Self::constant(1.0);
        }
        let abs_n = n.unsigned_abs() as usize;
        let mut result = *self;
        for _ in 1..abs_n {
            result = result * *self;
        }
        if n < 0 {
            Self::constant(1.0) / result
        } else {
            result
        }
    }

    /// Float power: exp(f * ln(self))
    fn powf(&self, f: f64) -> Self {
        (self.ln() * f).exp()
    }

    /// Jet power: exp(rhs * ln(self))
    fn pow(&self, rhs: Self) -> Self {
        (self.ln() * rhs).exp()
    }

    /// Square root using direct recurrence from $z^2 = a$:
    /// $z_0 = \sqrt{a_0}$,
    /// $z_n = \frac{1}{2\,z_0}\left(a_n - \sum_{k=1}^{n-1} z_k\, z_{n-k}\right)$
    fn sqrt(&self) -> Self {
        let a0 = self.coeff(0);
        let z0 = a0.sqrt();
        let inv_2z0 = 1.0 / (2.0 * z0);
        let mut z = Self::zero();
        z.set_coeff(0, z0);
        for n in 1..=N {
            let mut s = 0.0f64;
            for k in 1..n {
                s += z.coeff(k) * z.coeff(n - k);
            }
            z.set_coeff(n, inv_2z0 * (self.coeff(n) - s));
        }
        z
    }
}

// =============================================================================
// TrigOps
// =============================================================================

impl<const N: usize> TrigOps for Jet<N> {
    /// $\sin$ and $\cos$ computed together via coupled normalized recurrence:
    /// $s_0 = \sin(a_0)$, $c_0 = \cos(a_0)$,
    /// $s_n = \frac{1}{n}\sum_{k=1}^{n} k\, a_k\, c_{n-k}$,
    /// $c_n = -\frac{1}{n}\sum_{k=1}^{n} k\, a_k\, s_{n-k}$
    fn sin_cos(&self) -> (Self, Self) {
        let mut s = Self::zero();
        let mut c = Self::zero();
        s.set_coeff(0, self.coeff(0).sin());
        c.set_coeff(0, self.coeff(0).cos());
        for n in 1..=N {
            let mut ss = 0.0f64;
            let mut cs = 0.0f64;
            for k in 1..=n {
                let ka = (k as f64) * self.coeff(k);
                ss += ka * c.coeff(n - k);
                cs += ka * s.coeff(n - k);
            }
            s.set_coeff(n, ss / (n as f64));
            c.set_coeff(n, -cs / (n as f64));
        }
        (s, c)
    }

    fn sin(&self) -> Self {
        self.sin_cos().0
    }

    fn cos(&self) -> Self {
        self.sin_cos().1
    }

    fn tan(&self) -> Self {
        let (s, c) = self.sin_cos();
        s / c
    }

    /// $\sinh$ and $\cosh$ via coupled normalized recurrence (same as $\sin/\cos$ but no negative on $\cosh$):
    /// $s_n = \frac{1}{n}\sum_{k=1}^{n} k\, a_k\, c_{n-k}$,
    /// $c_n = \frac{1}{n}\sum_{k=1}^{n} k\, a_k\, s_{n-k}$
    fn sinh(&self) -> Self {
        self.sinh_cosh().0
    }

    fn cosh(&self) -> Self {
        self.sinh_cosh().1
    }

    fn tanh(&self) -> Self {
        let (s, c) = self.sinh_cosh();
        s / c
    }

    fn asin(&self) -> Self {
        // q = 1/sqrt(1 - a^2)
        let one = Self::constant(1.0);
        let q = (one - self.powi(2)).sqrt();
        let q_inv = one / q;
        self.integrate_derivative(self.coeff(0).asin(), &q_inv)
    }

    fn acos(&self) -> Self {
        // q = -1/sqrt(1 - a^2)
        let one = Self::constant(1.0);
        let q = (one - self.powi(2)).sqrt();
        let q_inv = -(one / q);
        self.integrate_derivative(self.coeff(0).acos(), &q_inv)
    }

    fn atan(&self) -> Self {
        // q = 1/(1 + a^2)
        let one = Self::constant(1.0);
        let q = one / (one + self.powi(2));
        self.integrate_derivative(self.coeff(0).atan(), &q)
    }

    fn asinh(&self) -> Self {
        // q = 1/sqrt(1 + a^2)
        let one = Self::constant(1.0);
        let q_inv = (one + self.powi(2)).sqrt();
        let q = one / q_inv;
        self.integrate_derivative(self.coeff(0).asinh(), &q)
    }

    fn acosh(&self) -> Self {
        // q = 1/sqrt(a^2 - 1)
        let one = Self::constant(1.0);
        let q_inv = (self.powi(2) - one).sqrt();
        let q = one / q_inv;
        self.integrate_derivative(self.coeff(0).acosh(), &q)
    }

    fn atanh(&self) -> Self {
        // q = 1/(1 - a^2)
        let one = Self::constant(1.0);
        let q = one / (one - self.powi(2));
        self.integrate_derivative(self.coeff(0).atanh(), &q)
    }
}

impl<const N: usize> Jet<N> {
    /// Compute $\sinh$ and $\cosh$ together via the coupled normalized recurrence.
    pub fn sinh_cosh(&self) -> (Self, Self) {
        let mut s = Self::zero();
        let mut c = Self::zero();
        s.set_coeff(0, self.coeff(0).sinh());
        c.set_coeff(0, self.coeff(0).cosh());
        for n in 1..=N {
            let mut ss = 0.0f64;
            let mut cs = 0.0f64;
            for k in 1..=n {
                let ka = (k as f64) * self.coeff(k);
                ss += ka * c.coeff(n - k);
                cs += ka * s.coeff(n - k);
            }
            s.set_coeff(n, ss / (n as f64));
            c.set_coeff(n, cs / (n as f64));  // NO negative for cosh
        }
        (s, c)
    }

    /// Integrate using derivative jet: used by inverse trig functions.
    /// Given $z'(a)$ encoded as a Jet `q`, compute $z$ coefficients by:
    /// $z_0 = z_0$,
    /// $z_n = \frac{1}{n}\sum_{k=1}^{n} k\, a_k\, q_{n-k}$
    fn integrate_derivative(&self, z0: f64, q: &Self) -> Self {
        let mut z = Self::zero();
        z.set_coeff(0, z0);
        for n in 1..=N {
            let mut s = 0.0f64;
            for k in 1..=n {
                s += (k as f64) * self.coeff(k) * q.coeff(n - k);
            }
            z.set_coeff(n, s / (n as f64));
        }
        z
    }
}

// =============================================================================
// ADFn — lift functions over Jet<2> to work at f64 or Jet level
// =============================================================================

/// Generic AD function wrapper.
///
/// Lifts a function `F: Fn(Jet<2>) -> Jet<2>` to operate at multiple levels:
/// - `call_stable(f64)` → `f64`: evaluate function value, first derivative, or second derivative
/// - `call_stable(Jet<2>)` → `Jet<2>`: pass through
///
/// For vector functions, also lifts `F: Fn(Vec<Jet<1>>) -> Vec<Jet<1>>`.
///
/// # Examples
/// ```
/// use peroxide::fuga::*;
///
/// let f_ad = ADFn::new(|x: Jet<2>| x.powi(2));
///
/// // Value: f(3) = 9
/// assert_eq!(f_ad.call_stable(3.0), 9.0);
///
/// // Gradient: f'(3) = 6  (2*3)
/// let df = f_ad.grad();
/// assert_eq!(df.call_stable(3.0), 6.0);
///
/// // Hessian: f''(3) = 2
/// let ddf = df.grad();
/// assert_eq!(ddf.call_stable(3.0), 2.0);
/// ```
pub struct ADFn<F> {
    f: Box<F>,
    grad_level: usize,
}

impl<F: Clone> ADFn<F> {
    /// Create a new `ADFn` wrapping function `f` at gradient level 0 (function evaluation).
    pub fn new(f: F) -> Self {
        Self {
            f: Box::new(f),
            grad_level: 0,
        }
    }

    /// Produce the gradient version of this function (increments grad_level by 1).
    /// Panics if grad_level >= 2.
    pub fn grad(&self) -> Self {
        assert!(self.grad_level < 2, "Higher order AD is not allowed");
        ADFn {
            f: self.f.clone(),
            grad_level: self.grad_level + 1,
        }
    }
}

/// Scalar version: F works with `Jet<2>`, target is `f64`.
impl<F: Fn(Jet<2>) -> Jet<2>> StableFn<f64> for ADFn<F> {
    type Output = f64;

    fn call_stable(&self, target: f64) -> f64 {
        match self.grad_level {
            0 => (self.f)(Jet::<2>::constant(target)).value(),
            1 => (self.f)(Jet::<2>::new(target, [1.0, 0.0])).dx(),
            2 => (self.f)(Jet::<2>::new(target, [1.0, 0.0])).ddx(),
            _ => unreachable!("grad_level > 2 is not allowed"),
        }
    }
}

/// Scalar version: F works with `Jet<2>`, target is `Jet<2>`.
impl<F: Fn(Jet<2>) -> Jet<2>> StableFn<Jet<2>> for ADFn<F> {
    type Output = Jet<2>;

    fn call_stable(&self, target: Jet<2>) -> Jet<2> {
        (self.f)(target)
    }
}

/// Vector version: F works with `Vec<Jet<1>>`, target is `Vec<f64>`.
impl<F: Fn(Vec<Jet<1>>) -> Vec<Jet<1>>> StableFn<Vec<f64>> for ADFn<F> {
    type Output = Vec<f64>;

    fn call_stable(&self, target: Vec<f64>) -> Vec<f64> {
        (self.f)(target.into_iter().map(Jet::<1>::constant).collect())
            .into_iter()
            .map(|j| j.value())
            .collect()
    }
}

/// Vector version: F works with `Vec<Jet<1>>`, target is `Vec<Jet<1>>`.
impl<F: Fn(Vec<Jet<1>>) -> Vec<Jet<1>>> StableFn<Vec<Jet<1>>> for ADFn<F> {
    type Output = Vec<Jet<1>>;

    fn call_stable(&self, target: Vec<Jet<1>>) -> Vec<Jet<1>> {
        (self.f)(target)
    }
}

/// Vector version: F works with `&Vec<Jet<1>>`, target is `&Vec<f64>`.
impl<'a, F: Fn(&Vec<Jet<1>>) -> Vec<Jet<1>>> StableFn<&'a Vec<f64>> for ADFn<F> {
    type Output = Vec<f64>;

    fn call_stable(&self, target: &'a Vec<f64>) -> Vec<f64> {
        let jet_target: Vec<Jet<1>> = target.iter().map(|&x| Jet::<1>::constant(x)).collect();
        (self.f)(&jet_target)
            .into_iter()
            .map(|j| j.value())
            .collect()
    }
}

/// Vector version: F works with `&Vec<Jet<1>>`, target is `&Vec<Jet<1>>`.
impl<'a, F: Fn(&Vec<Jet<1>>) -> Vec<Jet<1>>> StableFn<&'a Vec<Jet<1>>> for ADFn<F> {
    type Output = Vec<Jet<1>>;

    fn call_stable(&self, target: &'a Vec<Jet<1>>) -> Vec<Jet<1>> {
        (self.f)(target)
    }
}

// =============================================================================
// JetVec trait (replaces ADVec)
// =============================================================================

/// Trait for converting between `Vec<f64>` and `Vec<Jet<1>>`.
pub trait JetVec {
    fn to_jet_vec(&self) -> Vec<Jet<1>>;
    fn to_f64_vec(&self) -> Vec<f64>;
}

impl JetVec for Vec<f64> {
    fn to_jet_vec(&self) -> Vec<Jet<1>> {
        self.iter().map(|&x| Jet::<1>::constant(x)).collect()
    }

    fn to_f64_vec(&self) -> Vec<f64> {
        self.clone()
    }
}

impl JetVec for Vec<Jet<1>> {
    fn to_jet_vec(&self) -> Vec<Jet<1>> {
        self.clone()
    }

    fn to_f64_vec(&self) -> Vec<f64> {
        self.iter().map(|j| j.value()).collect()
    }
}

// =============================================================================
// FPVector, Vector, VecOps for Vec<Jet<1>>
// =============================================================================

impl FPVector for Vec<Jet<1>> {
    type Scalar = Jet<1>;

    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> Self::Scalar,
    {
        self.iter().map(|&x| f(x)).collect()
    }

    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
        T: Into<Self::Scalar>,
    {
        self.iter().fold(init.into(), |acc, &x| f(acc, x))
    }

    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
    {
        self.iter()
            .zip(other.iter())
            .map(|(&x, &y)| f(x, y))
            .collect()
    }

    fn filter<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> bool,
    {
        self.iter().filter(|&&x| f(x)).cloned().collect()
    }

    fn take(&self, n: usize) -> Self {
        self.iter().take(n).cloned().collect()
    }

    fn skip(&self, n: usize) -> Self {
        self.iter().skip(n).cloned().collect()
    }

    fn sum(&self) -> Self::Scalar {
        if self.is_empty() {
            return Jet::<1>::constant(0.0);
        }
        let s = self[0];
        self.reduce(s, |x, y| x + y)
    }

    fn prod(&self) -> Self::Scalar {
        if self.is_empty() {
            return Jet::<1>::constant(1.0);
        }
        let s = self[0];
        self.reduce(s, |x, y| x * y)
    }
}

impl Vector for Vec<Jet<1>> {
    type Scalar = Jet<1>;

    fn add_vec(&self, rhs: &Self) -> Self {
        self.add_v(rhs)
    }

    fn sub_vec(&self, rhs: &Self) -> Self {
        self.sub_v(rhs)
    }

    fn mul_scalar(&self, rhs: Self::Scalar) -> Self {
        self.mul_s(rhs)
    }
}

impl VecOps for Vec<Jet<1>> {}

// =============================================================================
// Backward compatibility: AD type aliases and constructors
// =============================================================================
// Keep the old public API so existing code that uses AD, AD0, AD1, AD2, ADVec, ADFn
// continues to compile. AD is now an alias for Jet<2> (the highest order used in ADFn).
// AD0/AD1/AD2 are re-exported constructor functions.

/// Backward compatibility alias: `AD` is now `Jet<2>`.
///
/// For new code, prefer `Dual = Jet<1>` or `HyperDual = Jet<2>` directly.
pub type AD = Jet<2>;

/// Backward compatibility constructor: `AD0(x)` creates a zero-derivative `Jet<2>` constant.
#[inline]
#[allow(non_snake_case)]
pub fn AD0(x: f64) -> Jet<2> {
    Jet::<2>::constant(x)
}

/// Backward compatibility constructor: `AD1(x, dx)` creates a `Jet<2>` with given first derivative.
#[inline]
#[allow(non_snake_case)]
pub fn AD1(x: f64, dx: f64) -> Jet<2> {
    Jet::<2>::new(x, [dx, 0.0])
}

/// Backward compatibility constructor: `AD2(x, dx, ddx)` creates a `Jet<2>` with given derivatives.
#[inline]
#[allow(non_snake_case)]
pub fn AD2(x: f64, dx: f64, ddx: f64) -> Jet<2> {
    Jet::<2>::new(x, [dx, ddx / 2.0])
}

/// Backward compatibility trait: provides `to_ad_vec` and `to_f64_vec` on vector types.
///
/// Extends `JetVec` with the `to_ad_vec` method for converting to `Vec<AD>` (= `Vec<Jet<2>>`).
pub trait ADVec: JetVec {
    fn to_ad_vec(&self) -> Vec<AD>;

    /// Convert to a `Vec<f64>` by extracting the value of each jet.
    /// (Delegates to `JetVec::to_f64_vec` — provided here so the trait is self-contained.)
    fn to_f64_vec_compat(&self) -> Vec<f64> {
        self.to_f64_vec()
    }
}

impl ADVec for Vec<f64> {
    fn to_ad_vec(&self) -> Vec<AD> {
        self.iter().map(|&x| Jet::<2>::constant(x)).collect()
    }
}

impl ADVec for Vec<AD> {
    fn to_ad_vec(&self) -> Vec<AD> {
        self.clone()
    }
}

impl JetVec for Vec<AD> {
    fn to_jet_vec(&self) -> Vec<Jet<1>> {
        self.iter()
            .map(|j| Jet::<1>::new(j.value(), [j.dx()]))
            .collect()
    }

    fn to_f64_vec(&self) -> Vec<f64> {
        self.iter().map(|j| j.value()).collect()
    }
}

// FPVector, Vector, VecOps for Vec<AD> (= Vec<Jet<2>>)
impl FPVector for Vec<AD> {
    type Scalar = AD;

    fn fmap<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> Self::Scalar,
    {
        self.iter().map(|&x| f(x)).collect()
    }

    fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
        T: Into<Self::Scalar>,
    {
        self.iter().fold(init.into(), |acc, &x| f(acc, x))
    }

    fn zip_with<F>(&self, f: F, other: &Self) -> Self
    where
        F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
    {
        self.iter()
            .zip(other.iter())
            .map(|(&x, &y)| f(x, y))
            .collect()
    }

    fn filter<F>(&self, f: F) -> Self
    where
        F: Fn(Self::Scalar) -> bool,
    {
        self.iter().filter(|&&x| f(x)).cloned().collect()
    }

    fn take(&self, n: usize) -> Self {
        self.iter().take(n).cloned().collect()
    }

    fn skip(&self, n: usize) -> Self {
        self.iter().skip(n).cloned().collect()
    }

    fn sum(&self) -> Self::Scalar {
        if self.is_empty() {
            return Jet::<2>::constant(0.0);
        }
        let s = self[0];
        self.reduce(s, |x, y| x + y)
    }

    fn prod(&self) -> Self::Scalar {
        if self.is_empty() {
            return Jet::<2>::constant(1.0);
        }
        let s = self[0];
        self.reduce(s, |x, y| x * y)
    }
}

impl Vector for Vec<AD> {
    type Scalar = AD;

    fn add_vec(&self, rhs: &Self) -> Self {
        self.add_v(rhs)
    }

    fn sub_vec(&self, rhs: &Self) -> Self {
        self.sub_v(rhs)
    }

    fn mul_scalar(&self, rhs: Self::Scalar) -> Self {
        self.mul_s(rhs)
    }
}

impl VecOps for Vec<AD> {}
