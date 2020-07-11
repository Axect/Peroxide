//! Easy to print any structures

use rand::distributions::uniform::SampleUniform;
use crate::statistics::dist::*;
use std::fmt::Debug;
use crate::structure::{
    matrix::Matrix,
    polynomial::Polynomial,
    multinomial::Multinomial,
    dual::{Dual, Dualist},
    hyper_dual::HyperDual,
    ad::AD,
};
use crate::traits::num::Number;
#[allow(unused_imports)]
#[cfg(feature = "dataframe")]
use crate::structure::dataframe::*;
use smallvec::Array;

pub trait Printable {
    fn print(&self);
}

impl Printable for f64 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for f32 {
    fn print(&self) {
        println!("{}", self);
    }
}
impl Printable for u64 {
    fn print(&self) {
        println!("{}", self);
    }
}
impl Printable for u32 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for usize {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for i64 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for i32 {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for Vec<usize> {
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

impl Printable for Vec<f64> {
    fn print(&self) {
        let mut result = String::new();
        result.push_str("[");
        for i in 0..self.len() {
            let st1 = format!("{:.4}", self[i]);
            let st2 = self[i].to_string();
            let mut st = st2.clone();

            if st1.len() < st2.len() {
                st = st1;
            }

            result.push_str(&st);
            if i == self.len() - 1 {
                break;
            }
            result.push_str(", ");
        }
        result.push_str("]");

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

impl Printable for Dual {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for Multinomial {
    fn print(&self) {
        println!("{}", self);
    }
}

impl Printable for Vec<Dual> {
    fn print(&self) {
        println!("value:");
        self.values().print();
        println!("slope:");
        self.slopes().print();
    }
}

impl Printable for HyperDual {
    fn print(&self) {
        println!("{}", self);
    }
}

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

impl Printable for Number {
    fn print(&self) {
        println!("{:?}", self)
    }
}

impl Printable for Vec<Number> {
    fn print(&self) {
        println!("{:?}", self);
    }
}

#[cfg(feature = "dataframe")]
impl Printable for DataFrame {
    fn print(&self) {
        println!("{}", self)
    }
}

impl<A: Array<Item=f64>> Printable for AD<A> {
    fn print(&self) {
        let mut result = String::new();
        result.push_str("AD [");
        for i in 0..self.len() {
            let st1 = format!("{:.4}", self[i]);
            let st2 = self[i].to_string();
            let mut st = st2.clone();

            if st1.len() < st2.len() {
                st = st1;
            }

            result.push_str(&st);
            if i == self.len() - 1 {
                break;
            }
            result.push_str(", ");
        }
        result.push_str("]");

        println!("{}", result);
    }
}
