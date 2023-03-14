//! Easy to print any structures

use crate::statistics::dist::*;
use crate::statistics::stat::ConfusionMatrix;
#[allow(unused_imports)]
use crate::structure::{
    ad::AD,
    matrix::Matrix,
    multinomial::Multinomial,
    polynomial::Polynomial,
    dataframe::{DataFrame, DTypeArray, Series, Scalar, DType},
};
use rand::distributions::uniform::SampleUniform;
use std::fmt::{Debug, LowerExp, UpperExp};

pub trait Printable {
    fn print(&self);
}

impl Printable for usize {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for u8 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for u16 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for u32 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for u64 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for isize {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for i8 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for i16 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for i32 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for i64 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for f32 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for f64 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for char {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for bool {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for String {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for Vec<usize> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<u8> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<u16> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<u32> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<u64> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<isize> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<i8> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<i16> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<i32> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<i64> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<char> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<&str> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for Vec<String> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

macro_rules! format_float_vec {
    ($self:expr) => {{
        let mut result = String::new();
        result.push_str("[");
        for i in 0 .. $self.len() {
            let st1 = $self[i].fmt_lower_exp(2);
            let st2 = $self[i].to_string();
            let st = if st1.len() < st2.len() {
                st1
            } else {
                st2
            };
            result.push_str(&st);
            if i == $self.len() - 1 {
                break;
            }
            result.push_str(", ");
        }
        result.push_str("]");
        result
    }};
}

impl Printable for Vec<f32> {
    fn print(&self) {
        let result = format_float_vec!(self);
        println!("{}", result);
    }
}

impl Printable for Vec<f64> {
    fn print(&self) {
        let result = format_float_vec!(self);
        println!("{}", result);
    }
}

impl Printable for &Vec<usize> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<u8> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<u16> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<u32> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<u64> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<isize> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<i8> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<i16> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<i32> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<i64> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<char> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<&str> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<String> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<bool> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl Printable for &Vec<f32> {
    fn print(&self) {
        let result = format_float_vec!(self);
        println!("{}", result);
    }
}

impl Printable for &Vec<f64> {
    fn print(&self) {
        let result = format_float_vec!(self);
        println!("{}", result);
    }
}

impl Printable for Matrix {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for Polynomial {
    fn print(&self) {
        println!("{}", self);
    }
}

//impl Printable for Dual {
//    fn print(&self) {
//        println!("{}", self);
//    }
//}

impl Printable for Multinomial {
    fn print(&self) {
        println!("{}", self);
    }
}

//impl Printable for Vec<Dual> {
//    fn print(&self) {
//        println!("value:");
//        self.values().print();
//        println!("slope:");
//        self.slopes().print();
//    }
//}
//
//impl Printable for HyperDual {
//    fn print(&self) {
//        println!("{}", self);
//    }
//}

impl<T: Debug + PartialOrd + SampleUniform + Copy + Into<f64>> Printable for OPDist<T> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

impl<T: Debug + PartialOrd + SampleUniform + Copy + Into<f64>> Printable for TPDist<T> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

//impl Printable for Number {
//    fn print(&self) {
//        println!("{:?}", self)
//    }
//}
//
//impl Printable for Vec<Number> {
//    fn print(&self) {
//        println!("{:?}", self);
//    }
//}

impl Printable for DType {
    fn print(&self) {
        println!("{}", self)
    }
}

impl Printable for DTypeArray {
    fn print(&self) {
        println!("{}", self)
    }
}

impl Printable for Scalar {
    fn print(&self) {
        println!("{}", self)
    }
}

impl Printable for Series {
    fn print(&self) {
        self.values.print();
    }
}

impl Printable for DataFrame {
    fn print(&self) {
        println!("{}", self)
    }
}

impl Printable for AD {
    fn print(&self) {
        println!("{}", self)
    }
}

impl Printable for ConfusionMatrix {
    fn print(&self) {
        println!("{}", self.to_matrix())
    }
}

/// Format float number into lower exponent notation with '+' sign
/// 
/// # Example
/// 
/// ```rust
/// use peroxide::fuga::*;
///
/// fn main() {
///     let x = 123.456;
///     assert_eq!(x.fmt_lower_exp(2), "1.23e+2");
/// }
/// ```
pub trait LowerExpWithPlus: LowerExp {
    fn fmt_lower_exp(&self, precision: usize) -> String {
        let mut s = format!("{:.p$e}", self, p=precision);
        let s_old = s.clone();
        let mut e = s.split_off(s.find('e').unwrap());
        if e.starts_with("e-") {
            s_old
        } else {
            e.insert(1, '+');
            format!("{}{}", s, e)
        }
    }
}

impl LowerExpWithPlus for f32 {}
impl LowerExpWithPlus for f64 {}

/// Format float number into upper exponent notation with '+' sign
/// 
/// # Example
/// 
/// ```rust
/// use peroxide::fuga::*;
///
/// fn main() {
///     let x = 123.456;
///     assert_eq!(x.fmt_upper_exp(2), "1.23E+2");
/// }
/// ```
pub trait UpperExpWithPlus: UpperExp {
    fn fmt_upper_exp(&self, precision: usize) -> String {
        let mut s = format!("{:.p$E}", self, p=precision);
        let s_old = s.clone();
        let mut e = s.split_off(s.find('E').unwrap());
        if e.starts_with("E-") {
            s_old
        } else {
            e.insert(1, '+');
            format!("{}{}", s, e)
        }
    }
}

impl UpperExpWithPlus for f32 {}
impl UpperExpWithPlus for f64 {}

//impl<A: Array<Item=f64>> Printable for AD<A> {
//    fn print(&self) {
//        let mut result = String::new();
//        result.push_str("AD [");
//        for i in 0..self.len() {
//            let st1 = format!("{:.4}", self[i]);
//            let st2 = self[i].to_string();
//            let mut st = st2.clone();
//
//            if st1.len() < st2.len() {
//                st = st1;
//            }
//
//            result.push_str(&st);
//            if i == self.len() - 1 {
//                break;
//            }
//            result.push_str(", ");
//        }
//        result.push_str("]");
//
//        println!("{}", result);
//    }
//}
