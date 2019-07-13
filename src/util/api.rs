use std::convert;
use structure::matrix::*;

pub trait MATLAB {
    fn new(s: &str) -> Self;
}

pub trait PYTHON {
    fn new<T>(v: Vec<Vec<T>>) -> Self
    where
        T: convert::Into<f64> + Copy;
}

pub trait R {
    fn new<T>(v: Vec<T>, x: usize, y: usize, shape: Shape) -> Self
    where
        T: convert::Into<f64>;
}

impl MATLAB for Matrix {
    fn new(s: &str) -> Self where {
        let str_rows_temp: Vec<&str> = s.split(';').collect();
        let str_rows = str_rows_temp
            .into_iter()
            .filter(|&t| t != "")
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
        Matrix {
            data: data,
            row: r,
            col: c,
            shape: Row,
        }
    }
}

impl PYTHON for Matrix {
    fn new<T>(v: Vec<Vec<T>>) -> Self
    where
        T: convert::Into<f64> + Copy,
    {
        let r = v.len();
        let c = v[0].len();
        let mut data = vec![0f64; r * c];
        for i in 0..r {
            for j in 0..c {
                let idx = i * c + j;
                data[idx] = v[i][j].into();
            }
        }
        Matrix {
            data: data,
            row: r,
            col: c,
            shape: Row,
        }
    }
}

impl R for Matrix {
    fn new<T>(v: Vec<T>, x: usize, y: usize, shape: Shape) -> Self
    where
        T: convert::Into<f64>,
    {
        Matrix {
            data: v.into_iter().map(|t| t.into()).collect::<Vec<f64>>(),
            row: x,
            col: y,
            shape: shape,
        }
    }
}
