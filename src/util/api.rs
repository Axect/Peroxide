//! Choose api - MATLAB, R, Python

use crate::structure::matrix::*;

pub trait MATLAB {
    fn new(s: &str) -> Self;
}

pub trait PYTHON {
    fn new<T>(v: Vec<Vec<T>>) -> Self
    where
        T: Into<f64> + Copy;
}

pub trait R {
    fn new<T>(v: Vec<T>, x: usize, y: usize, shape: Shape) -> Self
    where
        T: Into<f64>;
}

impl MATLAB for Matrix {
    fn new(s: &str) -> Self where {
        let str_rows_temp: Vec<&str> = s.split(';').collect();
        let str_rows = str_rows_temp
            .into_iter()
            .filter(|&t| !t.is_empty())
            .collect::<Vec<&str>>();
        let r = str_rows.len();
        let str_data = str_rows
            .into_iter()
            .map(|x| x.trim().split(' ').collect::<Vec<&str>>())
            .collect::<Vec<Vec<&str>>>();
        let c = str_data[0].len();
        let data = str_data
            .into_iter()
            .flat_map(|t| {
                t.into_iter()
                    .map(|x| x.parse::<f64>().unwrap())
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<f64>>();
        matrix(data, r, c, Row)
    }
}

impl PYTHON for Matrix {
    fn new<T>(v: Vec<Vec<T>>) -> Self
    where
        T: Into<f64> + Copy,
    {
        let r = v.len();
        let c = v[0].len();
        let mut data = vec![0f64; r * c];
        for (i, row) in v.iter().enumerate() {
            for (j, val) in row.iter().enumerate() {
                data[i * c + j] = (*val).into();
            }
        }
        matrix(data, r, c, Row)
    }
}

impl R for Matrix {
    fn new<T>(v: Vec<T>, x: usize, y: usize, shape: Shape) -> Self
    where
        T: Into<f64>,
    {
        matrix(
            v.into_iter().map(|t| t.into()).collect::<Vec<f64>>(),
            x,
            y,
            shape,
        )
    }
}
