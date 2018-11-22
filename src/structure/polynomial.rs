#[allow(unused_imports)]
use structure::matrix::*;
#[allow(unused_imports)]
use structure::vector::*;

use std::ops::{Neg, Add, Sub, Mul, Div};
use std::fmt;
use std::convert;
use std::cmp::{max, min};

// =============================================================================
// Polynomial Structure
// =============================================================================

/// Polynomial Structure
#[derive(Debug, Clone)]
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
            let temp = choose_shorter_string(
                format!("{}", value),
                format!("{:.4}", value),
            );
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
                let temp = choose_shorter_string(
                    format!("{}x", coef_1),
                    format!("{:.4}x", coef_1),
                );
                result.push_str(&temp);
            }

            if coef_0 > 0. {
                let temp = choose_shorter_string(
                    format!("+ {}", coef_0),
                    format!("+ {:.4}", coef_0),
                );
                result.push_str(&temp);
            } else if coef_0 < 0. {
                let temp = choose_shorter_string(
                    format!("- {}", coef_0.abs()),
                    format!("- {:.4}", coef_0.abs()),
                );
                result.push_str(&temp);
            }
            return write!(f, "{}", result);
        }

        for i in 0 .. l + 1 {
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
                },
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
                },
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
                },
                _ => {
                    let value = self.coef[i];
                    if value == 1. {
                        result.push_str(&format!(" + x^{}", l-i));
                    } else if value > 0. {
                        let temp = choose_shorter_string(
                            format!(" + {}x^{}", value, l-i),
                            format!(" + {:.4}x^{}", value, l-i),
                        );
                        result.push_str(&temp);
                    } else if value == -1. {
                        result.push_str(&format!(" - x^{}", l-i));
                    } else if value < 0. {
                        let temp = choose_shorter_string(
                            format!(" - {}x^{}", value.abs(), l-i),
                            format!(" - {:.4}x^{}", value.abs(), l-i),
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
    pub fn new(coef: Vector) -> Polynomial {
        Polynomial { coef }
    }

    /// Evaluate polynomial with value
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = poly(c!(1,3,2));
    /// assert_eq!(a.eval(1), 6_f64);
    /// ```
    pub fn eval<T>(&self, x: T) -> f64 where T: convert::Into<f64> + Copy {
        let l = self.coef.len() - 1;
        let mut s = 0f64;
        for i in 0 .. l + 1 {
            s += self.coef[i] * x.into().powf((l - i) as f64);
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
    type Output = Polynomial;

    fn neg(self) -> Self::Output {
        Polynomial::new(
            self.coef
                .clone()
                .into_iter()
                .map(|x| -x)
                .collect::<Vector>()
        )
    }
}

impl Add<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn add(self, other: Polynomial) -> Polynomial {
        let (l1, l2) = (self.coef.len(), other.coef.len());
        let l_max = max(l1, l2);
        let l_min = min(l1, l2);
        let v_max = choose_longer_vec(&self.coef, &other.coef);
        let v_min = choose_shorter_vec(&self.coef, &other.coef);
        let mut coef = vec![0f64; l_max];

        for i in 0 .. l_max {
            if i < l_max - l_min {
                coef[i] = v_max[i];
            } else {
                let j = i - (l_max - l_min);
                coef[i] = v_max[i] + v_min[j];
            }
        }
        Polynomial::new(coef)
    }
}

impl<T> Add<T> for Polynomial where T: convert::Into<f64> + Copy {
    type Output = Polynomial;
    fn add(self, other: T) -> Polynomial {
        Polynomial::new(self.coef.fmap(|x| x + other.into()))
    }
}

impl Sub<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn sub(self, other: Polynomial) -> Polynomial {
        self.add(other.neg())
    }
}

impl<T> Sub<T> for Polynomial where T: convert::Into<f64> + Copy {
    type Output = Polynomial;
    fn sub(self, other: T) -> Polynomial {
        Polynomial::new(self.coef.fmap(|x| x - other.into()))
    }
}

impl<T> Mul<T> for Polynomial 
    where T: convert::Into<f64> + Copy {
    type Output = Polynomial;
    fn mul(self, other: T) -> Polynomial {
        Polynomial::new(
            self.coef
                .into_iter()
                .map(|x| x * other.into())
                .collect::<Vector>()
        )
    }
}

impl Mul<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn mul(self, other: Polynomial) -> Polynomial {
        let (l1, l2) = (self.coef.len(), other.coef.len());
        let (n1, n2) = (l1 - 1, l2 - 1);
        let n = n1 + n2;
        let mut result = vec![0f64; n + 1];

        for i in 0 .. l1 {
            let fixed_val = self.coef[i];
            let fixed_exp = n1 - i;

            for j in 0 .. l2 {
                let target_val = other.coef[j];
                let target_exp = n2 - j;

                let result_val = fixed_val * target_val;
                let result_exp = fixed_exp + target_exp;

                result[n - result_exp] += result_val;
            }
        }

        Polynomial::new(result)
    }
}

impl<T> Div<T> for Polynomial where T: convert::Into<f64> + Copy {
    type Output = Polynomial;
    fn div(self, other: T) -> Polynomial {
        let val = other.into();
        assert_ne!(val, 0f64);

        Polynomial::new(self.coef.fmap(|x| x / val))
    }
}

impl Div<Polynomial> for Polynomial {
    type Output = (Polynomial, Polynomial);
    fn div(self, other: Polynomial) -> Self::Output {
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
            let mut temp_vec = vec![0f64; l-1];

            for i in 1 .. l {
                if i < l2 {
                    temp_vec[i-1] = temp.coef[i] - nom * other.coef[i];
                } else {
                    temp_vec[i-1] = temp.coef[i];
                }
            }

            temp = poly(temp_vec);
        }

        let rem = temp;

        (poly(quot_vec), rem)
    }
}

// =============================================================================
// Extra operations for Polynomial
// =============================================================================
pub trait ExtraOps {
    fn pow(&self, n: usize) -> Polynomial;
}

impl ExtraOps for Polynomial {
    fn pow(&self, n: usize) -> Polynomial {
        let mut result = self.clone();
        for _i in 0 .. n-1 {
            result = result * self.clone();
        }
        result
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

        for i in 0 .. l {
            result[i] = self.coef[i] * (l - i) as f64;
        }
        Polynomial::new(result)
    }

    fn integral(&self) -> Self {
        let l = self.coef.len();
        let mut result = vec![0f64; l + 1];

        for i in 0 .. l {
            result[i] = self.coef[i] / (l - i) as f64;
        }
        Polynomial::new(result)
    }
}


// =============================================================================
// Utils
// =============================================================================
fn choose_shorter_string(x1: String, x2: String) -> String {
    if x1.len() > x2.len() {
        x2
    } else {
        x1
    }
}

fn choose_shorter_vec(x1: &Vector, x2: &Vector) -> Vector {
    if x1.len() > x2.len() {
        x2.clone()
    } else {
        x1.clone()
    }
}

fn choose_longer_vec(x1: &Vector, x2: &Vector) -> Vector {
    if x1.len() <= x2.len() {
        x2.clone()
    } else {
        x1.clone()
    }
}
