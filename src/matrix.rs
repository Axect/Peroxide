use std::convert;
use std::fmt;
pub use self::Shape::{Row, Col};

#[derive(Debug)]
pub enum Shape {
    Row,
    Col,
}

#[derive(Debug)]
pub struct Matrix {
    data: Vec<f64>,
    row: usize,
    col: usize,
    shape: Shape,
}

trait Generic<T: convert::Into<f64>> {
    fn new(v: Vec<T>, x:usize, y:usize, shape: Shape) -> Matrix;
}

impl<T> Generic<T> for Matrix where T: convert::Into<f64> {
    fn new(v: Vec<T>, x: usize, y: usize, shape: Shape) -> Matrix {
        Matrix {
            data: v.into_iter().map(|x| x.into()).collect(),
            row: x,
            col: y,
            shape: shape,
        }
    }
}

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

#[allow(dead_code)]
impl Matrix {
    fn change_shape(&self) -> Matrix {
        let r = self.row;
        let c = self.col;
        assert_eq!(r*c, self.data.len());
        let l = r * c - 1;
        let mut data: Vec<f64> = self.data.clone();
        let ref_data: Vec<f64> = self.data.clone();

        match self.shape {
            Row => {
                for i in 0 .. l {
                    let s = (i * c) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                Matrix {
                    data: data.clone(),
                    row: r,
                    col: c,
                    shape: Col,
                }
            },
            Col => {
                for i in 0 .. l {
                    let s = (i * r) % l;
                    data[i] = ref_data[s];
                }
                data[l] = ref_data[l];
                Matrix {
                    data: data.clone(),
                    row: r,
                    col: c,
                    shape: Row,
                }
            },
        }
    }

    fn spread(&self) -> String {
        assert_eq!(self.row * self.col, self.data.len());
        let _rows = self.row;
        let cols = self.col;
        
        let mut result = String::new();
        result += &tab("");

        // Make header
        for i in 0 .. cols {
            let s = format!("c[{}]", i);
            result += &tab(&s);
            result += &s;
        }
        result += "\n";

        match self.shape {
            Row => {
                let data = self.data.clone();
                let temp: Vec<String> = data.into_iter().map(|x| x.to_string()).collect();
                let ts: Vec<String> = temp.clone().into_iter().take(cols).collect();
                let mut ss: Vec<String> = temp.into_iter().skip(cols).collect();
                let mut n: usize = 0;
                let s = format!("r[{}]", n);
                result += &s;
                result += &tab(&s);
                for txt in ts.into_iter() {
                    result += &tab(&txt);
                    result += &txt;
                }
                while ss.len() >= cols {
                    result += "\n";
                    let ts: Vec<String> = ss.clone().into_iter().take(cols).collect();
                    n += 1;
                    let s = format!("r[{}]", n);
                    result += &s;
                    result += &tab(&s);
                    let tl = ts.len();
                    for i in 0 .. tl {
                        result += &tab(&ts[i]);
                        result += &ts[i];
                    }
                    ss = ss.into_iter().skip(cols).collect();
                }
            },
            Col => {
                let mat = self.change_shape();
                let data = mat.data.clone();
                let temp: Vec<String> = data.into_iter().map(|x| x.to_string()).collect();
                let ts: Vec<String> = temp.clone().into_iter().take(cols).collect();
                let mut ss: Vec<String> = temp.into_iter().skip(cols).collect();
                let mut n: usize = 0;
                let s = format!("r[{}]", n);
                result += &s;
                result += &tab(&s);
                for txt in ts.into_iter() {
                    result += &tab(&txt);
                    result += &txt;
                }
                while ss.len() >= cols {
                    result += "\n";
                    let ts: Vec<String> = ss.clone().into_iter().take(cols).collect();
                    n += 1;
                    let s = format!("r[{}]", n);
                    result += &s;
                    result += &tab(&s);
                    let tl = ts.len();
                    for i in 0 .. tl {
                        result += &tab(&ts[i]);
                        result += &ts[i];
                    }
                    ss = ss.into_iter().skip(cols).collect();
                }
            }
        }
        return result;
    }
}

#[allow(unused_comparisons)]
pub fn tab(s: &str) -> String {
    let l = s.len();
    let mut m: String = String::new();
    if (5 - l) >= 0 {
        for _i in 0 .. (5 - l) {
            m += " ";
        }
    } else {
        m += " ";
    }
    return m;
}
