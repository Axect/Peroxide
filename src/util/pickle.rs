extern crate serde;
extern crate serde_pickle;

use std::fs::File;
use std::io::Write;
use std::process::exit;
use structure::matrix::*;
use structure::vector::*;

/// Pickle trait
///
/// # Description
///
/// Use python pickle to export vector or matrix
pub trait Pickle {
    fn write_single_pickle(&self, path: &str) -> serde_pickle::Result<()>;
    fn write_pickle<W: Write>(&self, writer: &mut W) -> serde_pickle::Result<()>;
}

impl Pickle for Vector {
    fn write_single_pickle(&self, path: &str) -> serde_pickle::Result<()> {
        let mut writer: Box<dyn Write>;

        match File::create(path) {
            Ok(p) => writer = Box::new(p),
            Err(e) => {
                println!("{:?}", e);
                exit(1);
            }
        }

        serde_pickle::to_writer(&mut writer, &self, true)
    }

    fn write_pickle<W: Write>(&self, writer: &mut W) -> serde_pickle::Result<()> {
        serde_pickle::to_writer(writer, &self, true)
    }
}

impl Pickle for Matrix {
    fn write_single_pickle(&self, path: &str) -> serde_pickle::Result<()> {
        let mut writer: Box<dyn Write>;

        match File::create(path) {
            Ok(p) => writer = Box::new(p),
            Err(e) => {
                println!("{:?}", e);
                exit(1);
            }
        }

        let mut container: Vec<Vec<f64>> = Vec::new();;

        match self.shape {
            Row => {
                for i in 0..self.row {
                    container.push(self.row(i));
                }
            }
            Col => {
                for i in 0..self.col {
                    container.push(self.col(i));
                }
            }
        }

        serde_pickle::to_writer(&mut writer, &container, true)
    }

    fn write_pickle<W: Write>(&self, writer: &mut W) -> serde_pickle::Result<()> {
        let mut container: Vec<Vec<f64>> = Vec::new();;

        match self.shape {
            Row => {
                for i in 0..self.row {
                    container.push(self.row(i));
                }
            }
            Col => {
                for i in 0..self.col {
                    container.push(self.col(i));
                }
            }
        }

        serde_pickle::to_writer(writer, &container, true)
    }
}
