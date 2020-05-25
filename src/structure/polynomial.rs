use operation::extra_ops::PowOps;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
#[allow(unused_imports)]
use structure::matrix::*;
#[allow(unused_imports)]
use structure::vector::*;
use util::useful::*;

use std::cmp::{max, min};
use std::convert;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

// =============================================================================
// Polynomial Structure
// =============================================================================

/// Polynomial Structure
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polynomial {
    pub coef: Vector,
}

/// Polynomial Print
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = poly(c!(1,3,2));
/// a.print(); //x^2 + 3x + 2
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
    pub fn new(coef: Vector) -> Self {
        Self { coef }
    }

    /// Evaluate polynomial with value according to Horner's method
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = poly(c!(1,3,2));
    /// assert_eq!(a.eval(1), 6_f64);
    ///
    /// let b = poly(c!(1, 1, -2, -2));
    /// let x = 2_f64.sqrt();
    /// let horner_evaluation = b.eval(x);
    /// let naive_evaluation = x.powf(3.0) + x.powf(2.0) - 2.0*x - 2.0;
    /// assert_eq!(horner_evaluation, 0_f64);
    /// assert_ne!(naive_evaluation, horner_evaluation);
    /// ```
    pub fn eval<T>(&self, x: T) -> f64
    where
        T: convert::Into<f64> + Copy,
    {
        let l = self.coef.len() - 1;
        let mut s = self.coef[0];
        for i in 0..l {
            s = self.coef[i+1] + x.into()*s;
        }
        s
    }

    pub fn eval_vec(&self, v: Vector) -> Vector {
        v.fmap(|t| self.eval(t))
    }
}

/// Convenient to declare polynomial
pub fn poly(coef: Vector) -> Polynomial {
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
                .collect::<Vector>(),
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
                .collect::<Vector>(),
        )
    }
}

impl Mul<Polynomial> for Polynomial {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
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
        let mut quot_vec: Vector = Vec::new();
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
pub fn lagrange_polynomial(node_x: Vector, node_y: Vector) -> Polynomial {
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