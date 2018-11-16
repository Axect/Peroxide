#[allow(unused_imports)]
use structure::matrix::*;
#[allow(unused_imports)]
use structure::vector::*;
#[allow(unused_imports)]
use structure::polynomial::*;

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

impl Printable for Vector {
    fn print(&self) {
        let mut result = String::new();
        result.push_str("[");
        for i in 0 .. self.len() {
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