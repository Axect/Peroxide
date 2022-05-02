#[allow(unused_imports)]
use crate::structure::matrix::*;
#[allow(unused_imports)]
use crate::structure::vector::*;
use crate::util::useful::*;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::traits::{fp::FPVector, num::PowOps};
use std::cmp::{max, min};
use std::convert;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

// =============================================================================
// Polynomial Structure
// =============================================================================
/// Polynomial Structure
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polynomial {
    pub coef: Vec<f64>,
}

/// Polynomial Print
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let a = poly(c!(1,3,2));
///     a.print(); //x^2 + 3x + 2
/// }
/// ```
impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut result = String::new();
        let l = self.coef.len() - 1;

        if l == 0 {
            let value = self.coef[0];
            let temp = choose_shorter_string(format!("{}", value), format!("{:.4}", value));
            return write!(f, "{}", temp);
        }

        if l == 1 {
            let coef_1 = self.coef[0];
            let coef_0 = self.coef[1];

            if coef_1 == 1. {
                result.push_str("x");
            } else if coef_1 == -1. {
                result.push_str("-x");
            } else {
                let temp = choose_shorter_string(format!("{}x", coef_1), format!("{:.4}x", coef_1));
                result.push_str(&temp);
            }

            if coef_0 > 0. {
                let temp =
                    choose_shorter_string(format!(" + {}", coef_0), format!(" + {:.4}", coef_0));
                result.push_str(&temp);
            } else if coef_0 < 0. {
                let temp = choose_shorter_string(
                    format!(" - {}", coef_0.abs()),
                    format!(" - {:.4}", coef_0.abs()),
                );
                result.push_str(&temp);
            }
            return write!(f, "{}", result);
        }

        for i in 0..l + 1 {
            match i {
                0 => {
                    let value = self.coef[i];
                    if value == 1. {
                        result.push_str(&format!("x^{}", l));
                    } else if value == -1. {
                        result.push_str(&format!("-x^{}", l));
                    } else {
                        let temp = choose_shorter_string(
                            format!("{}x^{}", value, l),
                            format!("{:.4}x^{}", value, l),
                        );
                        result.push_str(&temp);
                    }
                }
                i if i == l => {
                    let value = self.coef[i];
                    if value > 0. {
                        let temp = choose_shorter_string(
                            format!(" + {}", value),
                            format!(" + {:.4}", value),
                        );
                        result.push_str(&temp);
                    } else if value < 0. {
                        let temp = choose_shorter_string(
                            format!(" - {}", value.abs()),
                            format!(" - {:.4}", value.abs()),
                        );
                        result.push_str(&temp);
                    }
                }
                i if i == l - 1 => {
                    let value = self.coef[i];
                    if value == 1. {
                        result.push_str(" + x");
                    } else if value > 0. {
                        let temp = choose_shorter_string(
                            format!(" + {}x", value),
                            format!(" + {:.4}x", value),
                        );
                        result.push_str(&temp);
                    } else if value == -1. {
                        result.push_str(" - x");
                    } else if value < 0. {
                        let temp = choose_shorter_string(
                            format!(" - {}x", value.abs()),
                            format!(" - {:.4}x", value.abs()),
                        );
                        result.push_str(&temp);
                    }
                }
                _ => {
                    let value = self.coef[i];
                    if value == 1. {
                        result.push_str(&format!(" + x^{}", l - i));
                    } else if value > 0. {
                        let temp = choose_shorter_string(
                            format!(" + {}x^{}", value, l - i),
                            format!(" + {:.4}x^{}", value, l - i),
                        );
                        result.push_str(&temp);
                    } else if value == -1. {
                        result.push_str(&format!(" - x^{}", l - i));
                    } else if value < 0. {
                        let temp = choose_shorter_string(
                            format!(" - {}x^{}", value.abs(), l - i),
                            format!(" - {:.4}x^{}", value.abs(), l - i),
                        );
                        result.push_str(&temp);
                    }
                }
            }
        }
        write!(f, "{}", result)
    }
}

impl Polynomial {
    /// Create Polynomial
    pub fn new(coef: Vec<f64>) -> Self {
        Self { coef }
    }

    /// Evaluate polynomial with value according to Horner's method
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = poly(c!(1,3,2));
    ///     assert_eq!(a.eval(1), 6_f64);
    ///
    ///     let b = poly(c!(1, 1, -2, -2));
    ///     let x = 2_f64.sqrt();
    ///     let horner_evaluation = b.eval(x);
    ///     let naive_evaluation = x.powf(3.0) + x.powf(2.0) - 2.0*x - 2.0;
    ///     assert_eq!(horner_evaluation, 0_f64);
    ///     assert_ne!(naive_evaluation, horner_evaluation);
    /// }
    /// ```
    pub fn eval<T>(&self, x: T) -> f64
    where
        T: convert::Into<f64> + Copy,
    {
        let x = x.into();
        let l = self.coef.len() - 1;
        let mut s = self.coef[0];
        for i in 0..l {
            s = self.coef[i + 1] + x * s;
        }
        s
    }

    pub fn eval_vec(&self, v: Vec<f64>) -> Vec<f64> {
        v.fmap(|t| self.eval(t))
    }

    /// Linear transformation of a polynomial by a given x according to Horner's method
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = poly(c!(1,3,2));
    ///     let translated = a.translate_x(2);
    ///
    ///     assert_eq!(translated.eval(3), 6_f64);
    /// }
    /// ```
    pub fn translate_x<X>(&self, x: X) -> Self
    where
        X: convert::Into<f64> + Copy,
    {
        let d = Self::new(vec![1f64, x.into()]);

        let mut coef = vec![0f64; self.coef.len()];

        let (mut p, ri) = self.honor_division(&d);
        coef[self.coef.len() - 1] = ri;

        for i in (0..(self.coef.len() - 1)).rev() {
            if p.coef.len() == 1 {
                coef[i] = p.coef[0];
                break;
            }

            let t = p.honor_division(&d);
            coef[i] = t.1;
            p = t.0;
        }

        Self::new(coef)
    }

    fn honor_division(&self, other: &Self) -> (Self, f64) {
        assert_eq!(other.coef.len(), 2usize);
        assert_eq!(other.coef[0], 1.0f64);

        let mut coef = vec![0f64; self.coef.len() - 1];
        coef[0] = self.coef[0];

        let d = other.coef[1];
        for i in 1..coef.len() {
            coef[i] = self.coef[i] - d * coef[i - 1];
        }

        let remainder = self.coef[self.coef.len() - 1] - d * coef[coef.len() - 1];
        (Self::new(coef), remainder)
    }
}

fn polynomial_fft_mul(poly1: &Polynomial, poly2: &Polynomial) -> Polynomial {
    use crate::numerical::fft;
    let target_length = Some(2 * usize::max(poly1.coef.len(), poly2.coef.len()));

    let prod_length = poly1.coef.len() + poly2.coef.len() - 1;

    let mut poly1_transform = fft::fft(&poly1.coef, target_length);
    let poly2_transform = fft::fft(&poly2.coef, target_length);

    for (val1, val2) in poly1_transform.iter_mut().zip(poly2_transform.iter()) {
        *val1 *= val2;
    }

    let mut res = fft::inverse_fft(&poly1_transform);

    res.resize(prod_length, 0.);

    poly(res)
}

/// Convenient to declare polynomial
pub fn poly(coef: Vec<f64>) -> Polynomial {
    Polynomial::new(coef)
}

// =============================================================================
// std::ops for Polynomial
// =============================================================================

impl Neg for Polynomial {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(
            self.coef
                .clone()
                .into_iter()
                .map(|x| -x)
                .collect::<Vec<f64>>(),
        )
    }
}

impl Add<Polynomial> for Polynomial {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let (l1, l2) = (self.coef.len(), other.coef.len());
        let l_max = max(l1, l2);
        let l_min = min(l1, l2);
        let v_max = choose_longer_vec(&self.coef, &other.coef);
        let v_min = choose_shorter_vec(&self.coef, &other.coef);
        let mut coef = vec![0f64; l_max];

        for i in 0..l_max {
            if i < l_max - l_min {
                coef[i] = v_max[i];
            } else {
                let j = i - (l_max - l_min);
                coef[i] = v_max[i] + v_min[j];
            }
        }
        Self::new(coef)
    }
}

impl<T> Add<T> for Polynomial
where
    T: convert::Into<f64> + Copy,
{
    type Output = Self;
    fn add(self, other: T) -> Self {
        Self::new(self.coef.fmap(|x| x + other.into()))
    }
}

impl Sub<Polynomial> for Polynomial {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self.add(other.neg())
    }
}

impl<T> Sub<T> for Polynomial
where
    T: convert::Into<f64> + Copy,
{
    type Output = Self;
    fn sub(self, other: T) -> Self {
        Self::new(self.coef.fmap(|x| x - other.into()))
    }
}

impl<T> Mul<T> for Polynomial
where
    T: convert::Into<f64> + Copy,
{
    type Output = Self;
    fn mul(self, other: T) -> Self {
        Self::new(
            self.coef
                .into_iter()
                .map(|x| x * other.into())
                .collect::<Vec<f64>>(),
        )
    }
}

const SWITCH_T0_FFT_FOR_MUL: usize = 1000;

impl Mul<Polynomial> for Polynomial {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        if usize::max(self.coef.len(), other.coef.len()) < SWITCH_T0_FFT_FOR_MUL {
            let (l1, l2) = (self.coef.len(), other.coef.len());
            let (n1, n2) = (l1 - 1, l2 - 1);
            let n = n1 + n2;
            let mut result = vec![0f64; n + 1];

            for i in 0..l1 {
                let fixed_val = self.coef[i];
                let fixed_exp = n1 - i;

                for j in 0..l2 {
                    let target_val = other.coef[j];
                    let target_exp = n2 - j;

                    let result_val = fixed_val * target_val;
                    let result_exp = fixed_exp + target_exp;

                    result[n - result_exp] += result_val;
                }
            }

            Self::new(result)
        } else {
            polynomial_fft_mul(&self, &other)
        }
    }
}

impl<T> Div<T> for Polynomial
where
    T: convert::Into<f64> + Copy,
{
    type Output = Self;
    fn div(self, other: T) -> Self {
        let val = other.into();
        assert_ne!(val, 0f64);

        Self::new(self.coef.fmap(|x| x / val))
    }
}

impl Div<Polynomial> for Polynomial {
    type Output = (Self, Self);
    fn div(self, other: Self) -> Self::Output {
        let l1 = self.coef.len();
        let l2 = other.coef.len();
        assert!(l1 >= l2);

        let mut temp = self.clone();
        let mut quot_vec: Vec<f64> = Vec::new();
        let denom = other.coef[0];

        while temp.coef.len() >= l2 {
            let l = temp.coef.len();
            let target = temp.coef[0];
            let nom = target / denom;
            quot_vec.push(nom);
            let mut temp_vec = vec![0f64; l - 1];

            for i in 1..l {
                if i < l2 {
                    temp_vec[i - 1] = temp.coef[i] - nom * other.coef[i];
                } else {
                    temp_vec[i - 1] = temp.coef[i];
                }
            }

            temp = poly(temp_vec);
        }

        let rem = temp;

        (poly(quot_vec), rem)
    }
}

impl Mul<Polynomial> for usize {
    type Output = Polynomial;

    fn mul(self, rhs: Polynomial) -> Self::Output {
        rhs.mul(self as f64)
    }
}

impl Mul<Polynomial> for i32 {
    type Output = Polynomial;

    fn mul(self, rhs: Polynomial) -> Self::Output {
        rhs.mul(self as f64)
    }
}

impl Mul<Polynomial> for i64 {
    type Output = Polynomial;

    fn mul(self, rhs: Polynomial) -> Self::Output {
        rhs.mul(self as f64)
    }
}

impl Mul<Polynomial> for f32 {
    type Output = Polynomial;

    fn mul(self, rhs: Polynomial) -> Self::Output {
        rhs.mul(self as f64)
    }
}

impl Mul<Polynomial> for f64 {
    type Output = Polynomial;

    fn mul(self, rhs: Polynomial) -> Self::Output {
        rhs.mul(self as f64)
    }
}

// =============================================================================
// Extra operations for Polynomial
// =============================================================================
impl PowOps for Polynomial {
    fn powi(&self, n: i32) -> Self {
        let mut result = self.clone();
        for _i in 0..n - 1 {
            result = result * self.clone();
        }
        result
    }

    fn powf(&self, _f: f64) -> Self {
        unimplemented!()
    }

    fn pow(&self, _f: Self) -> Self {
        unimplemented!()
    }

    fn sqrt(&self) -> Self {
        unimplemented!()
    }
}

// =============================================================================
// Calculus for Polynomial
// =============================================================================
pub trait Calculus {
    fn diff(&self) -> Self;
    fn integral(&self) -> Self;
}

impl Calculus for Polynomial {
    fn diff(&self) -> Self {
        let l = self.coef.len() - 1;
        let mut result = vec![0f64; l];

        for i in 0..l {
            result[i] = self.coef[i] * (l - i) as f64;
        }
        Self::new(result)
    }

    fn integral(&self) -> Self {
        let l = self.coef.len();
        let mut result = vec![0f64; l + 1];

        for i in 0..l {
            result[i] = self.coef[i] / (l - i) as f64;
        }
        Self::new(result)
    }
}

// =============================================================================
// Useful Polynomial
// =============================================================================
/// Lagrange Polynomial
pub fn lagrange_polynomial(node_x: Vec<f64>, node_y: Vec<f64>) -> Polynomial {
    assert_eq!(node_x.len(), node_y.len());
    let l = node_x.len();
    let mut result = Polynomial::new(vec![0f64; l]);

    for i in 0..l {
        let fixed_val = node_x[i];
        let prod = node_y[i];
        let mut id = poly(vec![1f64]);

        for j in 0..l {
            if j == i {
                continue;
            } else {
                let target_val = node_x[j];
                let denom = fixed_val - target_val;
                id = id * (poly(vec![1f64, -target_val]) / denom);
            }
        }
        result = result + (id * prod);
    }
    result
}

/// Legendre Polynomial
///
/// # Description
/// : Generate `n`-th order of Legendre polynomial
pub fn legendre_polynomial(n: usize) -> Polynomial {
    match n {
        0 => poly(vec![1f64]),       // 1
        1 => poly(vec![1f64, 0f64]), // x
        2 => poly(vec![1.5, 0f64, -0.5]),
        3 => poly(vec![2.5, 0f64, -1.5, 0f64]),
        _ => {
            let k = n - 1;
            let k_f64 = k as f64;
            ((2f64 * k_f64 + 1f64) * poly(vec![1f64, 0f64]) * legendre_polynomial(k)
                - k_f64 * legendre_polynomial(k - 1))
                / (k_f64 + 1f64)
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecialKind {
    First,
    Second
}

/// Chebyshev Polynomial
pub fn chebyshev_polynomial(n: usize, kind: SpecialKind) -> Polynomial {
    let mut prev = Polynomial::new(vec![1f64]);
    let mut curr = match kind {
        SpecialKind::First => Polynomial::new(vec![1f64, 0f64]),
        SpecialKind::Second => Polynomial::new(vec![2f64, 0f64]),
    };

    match n {
        0 => prev,
        1 => curr,
        _ => {
            for _i in 1 .. n {
                std::mem::swap(&mut prev, &mut curr);
                curr = poly(vec![2f64, 0f64]) * prev.clone() - curr;
            }
            curr
        }
    }
}

#[cfg(test)]
mod tests {
    use std::task::Poll;

    use super::*;
    use float_cmp::approx_eq;

    #[test]
    fn test_honor_division() {
        let a = Polynomial::new(vec![1f64, -4f64, 4f64, 3f64, -8f64, 4f64]);
        let b = Polynomial::new(vec![1f64, -2f64]);

        let (c, remainder) = a.honor_division(&b);
        assert_eq!(c.coef, vec![1f64, -2f64, 0f64, 3f64, -2f64]);
        assert_eq!(remainder, 0f64);
    }

    #[test]
    fn test_translate_x() {
        let a = Polynomial::new(vec![1f64, -4f64, 4f64, 3f64, -8f64, 4f64]);
        let b = a.translate_x(-6);

        for i in -10..10 {
            assert_eq!(a.eval(i), b.eval(i - 6));
        }
    }

    #[test]
    fn test_low_degree_mul() {
        let a = Polynomial::new(vec![1., 1.]);
        let b = Polynomial::new(vec![-1., 1.]);

        let target = Polynomial::new(vec![-1., 0., 1.]);
        let prod = a * b;

        for (prod_coef, target_coef) in target.coef.iter().zip(prod.coef.iter()) {
            assert!(approx_eq!(f64, *prod_coef, *target_coef));
        }
    }

    #[test]
    fn test_constant_poly_mul() {
        let a = Polynomial::new(vec![3.]);
        let b = Polynomial::new(vec![1., 2., 3., 4., 5.]);

        let prod = a * b.clone();

        for (prod_coef, orig_coef) in prod.coef.iter().zip(b.coef.iter()) {
            assert!(approx_eq!(f64, *prod_coef, 3. * *orig_coef));
        }
    }

    #[test]
    fn test_large_poly_mul() {
        let a = Polynomial::new(
            (1..=SWITCH_T0_FFT_FOR_MUL)
                .into_iter()
                .map(|i| i as f64)
                .collect(),
        );
        let b = Polynomial::new(vec![1.; SWITCH_T0_FFT_FOR_MUL]);

        let target_coefs = (1..2 * (SWITCH_T0_FFT_FOR_MUL - 1))
            .into_iter()
            .map(|i| {
                if i <= 1000 {
                    (i * (i + 1) / 2) as f64
                } else {
                    (1000 * 1001 / 2 - (i - 1000) * (i - 999) / 2) as f64
                }
            })
            .collect::<Vec<f64>>();

        for (i, (coef, target_coef)) in (a * b).coef.iter().zip(target_coefs.iter()).enumerate() {
            println!("{} {} {}", i, coef, target_coef);
            assert!(approx_eq!(
                f64,
                *coef,
                *target_coef,
                epsilon = (10.).powi(-8)
            ));
        }
    }
}
