#[allow(unused_imports)]
use matrix::*;
#[allow(unused_imports)]
use vector::*;

use std::ops::{Neg, Add, Sub, Mul};
use std::fmt;
use std::convert;
use std::cmp::{max, min};

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coef: Vector,
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut result = String::new();
        let l = self.coef.len() - 1;

        for i in 0 .. l + 1 {
            match i {
                0 => {
                    let value = self.coef[i];
                    if value == 1. {
                        result.push_str(&format!("x^{}", l));
                    } else if value == -1. {
                        result.push_str(&format!("-x^{}", l));
                    } else {
                        let temp = choose_shorter(
                        format!("{}x^{}", value, l),
                        format!("{:.4}x^{}", value, l),
                        );
                        result.push_str(&temp);
                    }
                },
                i if i == l => {
                    let value = self.coef[i];
                    if value > 0. {
                        let temp = choose_shorter(
                            format!(" + {}", value),
                            format!(" + {:.4}", value),
                        );
                        result.push_str(&temp);
                    } else if value < 0. {
                        let temp = choose_shorter(
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
                        let temp = choose_shorter(
                            format!(" + {}x", value),
                            format!(" + {:.4}x", value),
                        );
                        result.push_str(&temp);
                    } else if value == -1. {
                        result.push_str(" - x");
                    } else if value < 0. {
                        let temp = choose_shorter(
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
                        let temp = choose_shorter(
                            format!(" + {}x^{}", value, l-i),
                            format!(" + {:.4}x^{}", value, l-i),
                        );
                        result.push_str(&temp);
                    } else if value == -1. {
                        result.push_str(&format!(" - x^{}", l-i));
                    } else if value < 0. {
                        let temp = choose_shorter(
                            format!(" - {}x^{}", value.abs(), l-i),
                            format!(" - {:.4}x{}", value.abs(), l-i),
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
    pub fn new(coef: Vector) -> Polynomial {
        Polynomial { coef: coef }
    }

    pub fn eval<T>(&self, x: T) -> f64 where T: convert::Into<f64> + Copy {
        let l = self.coef.len() - 1;
        let mut s = 0f64;
        for i in 0 .. l + 1 {
            s += self.coef[i] * x.into().powf((l - i) as f64);
        }
        s
    } 
}

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

impl Sub<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn sub(self, other: Polynomial) -> Polynomial {
        self.add(other.neg())
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

pub struct Lagrange {
    x: Vector,
    y: Vector,
}

impl Lagrange {
    pub fn new(x: Vector, y: Vector) -> Lagrange {
        Lagrange { x: x, y: y }
    }
}

fn choose_shorter(x1: String, x2: String) -> String {
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
    if x1.len() < x2.len() {
        x2.clone()
    } else {
        x1.clone()
    }
}