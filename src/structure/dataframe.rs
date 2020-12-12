//use self::csv::{ReaderBuilder, WriterBuilder};
//use crate::structure::matrix::{matrix, Matrix, Shape::*};
//use crate::traits::fp::FPMatrix;
//use crate::util::useful::tab;
//use indexmap::{map::Keys, IndexMap};
//use json::JsonValue;
//use std::cmp::{max, min};
//use std::collections::HashMap;
//use std::error::Error;
use std::fmt;
use std::ops::{Index, IndexMut};
use std::cmp::{max, min};
use crate::util::useful::tab;
use DType::{
    USIZE,U8,U16,U32,U64,
    ISIZE,I8,I16,I32,I64,
    F32,F64,Bool,Char,Str
};
//use std::{fmt, fmt::Debug};

// =============================================================================
// Enums
// =============================================================================

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum DType {
    USIZE,
    U8,
    U16,
    U32,
    U64,
    ISIZE,
    I8,
    I16,
    I32,
    I64,
    F32,
    F64,
    Bool,
    Str,
    Char,
}

#[derive(Debug, Clone)]
pub enum DTypeArray {
    USIZE(Vec<usize>),
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
    U64(Vec<u64>),
    ISIZE(Vec<isize>),
    I8(Vec<i8>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    I64(Vec<i64>),
    F32(Vec<f32>),
    F64(Vec<f64>),
    Bool(Vec<bool>),
    Str(Vec<String>),
    Char(Vec<char>),
}

#[derive(Debug, Clone)]
pub enum DTypeValue {
    USIZE(usize),
    U8(u8),
    U16(u16),
    U32(u32),
    U64(u64),
    ISIZE(isize),
    I8(i8),
    I16(i16),
    I32(i32),
    I64(i64),
    F32(f32),
    F64(f64),
    Bool(bool),
    Str(String),
    Char(char),
}

// =============================================================================
// Structs
// =============================================================================

#[derive(Debug, Clone)]
pub struct DataFrame {
    pub data: Vec<Series>,
    pub ics: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct Series {
    pub values: DTypeArray,
    pub dtype: DType,
}

#[derive(Debug, Clone)]
pub struct Scalar {
    pub value: DTypeValue,
    pub dtype: DType,
}

pub struct IndexSeries {
    pub name: String,
    pub ics: Vec<String>,
}

// =============================================================================
// Traits
// =============================================================================
pub trait TypedScalar<T> {
    fn new(s: T) -> Self where Self: Sized;
    fn unwrap(self) -> T;
}

pub trait TypedVector<T> {
    fn new(v: Vec<T>) -> Self;
    fn to_vec(self) -> Vec<T>;
    fn as_slice(&self) -> &[T];
    fn at_raw(&self, i: usize) -> T;
}


// =============================================================================
// Macros
// =============================================================================
macro_rules! impl_typed_scalar {
    ($type:ty, $dtype:ident) => {
        impl TypedScalar<$type> for Scalar {
            fn new(s: $type) -> Self {
                Self {
                    value: DTypeValue::$dtype(s),
                    dtype: DType::$dtype,
                }
            }

            fn unwrap(self) -> $type {
                match self.value {
                    DTypeValue::$dtype(s) => s,
                    _ => panic!("Can't unwrap {:?} value", $dtype),
                }
            }
        }
    }
}

macro_rules! impl_typed_vector {
    ($type:ty, $dtype:ident) => {
        impl TypedVector<$type> for Series {
            fn new(v: Vec<$type>) -> Self {
                Self {
                    values: DTypeArray::$dtype(v),
                    dtype: DType::$dtype,
                }
            }

            fn to_vec(self) -> Vec<$type> {
                match self.values {
                    DTypeArray::$dtype(v) => v,
                    _ => panic!("Can't convert to {:?} vector", $dtype),
                }
            }

            fn as_slice(&self) -> &[$type] {
                match &self.values {
                    DTypeArray::$dtype(v) => v,
                    _ => panic!("Can't convert to {:?} vector", $dtype),
                }
            }

            fn at_raw(&self, i: usize) -> $type {
                let v: &[$type] = self.as_slice();
                v[i].clone()
            }
        }
    }
}

macro_rules! dtype_case {
    ($type:ty, $value:expr, $wrapper: expr) => {{
        let x: $type = $value;
        $wrapper(x)
    }};
}

macro_rules! dtype_match {
    ($dtype:expr, $value:expr, $wrapper:expr) => {{
        match $dtype {
            USIZE => dtype_case!(usize, $value, $wrapper),
            U8 => dtype_case!(u8, $value, $wrapper),
            U16 => dtype_case!(u16, $value, $wrapper),
            U32 => dtype_case!(u32, $value, $wrapper),
            U64 => dtype_case!(u64, $value, $wrapper),
            ISIZE => dtype_case!(isize, $value, $wrapper),
            I8 => dtype_case!(i8, $value, $wrapper),
            I16 => dtype_case!(i16, $value, $wrapper),
            I32 => dtype_case!(i32, $value, $wrapper),
            I64 => dtype_case!(i64, $value, $wrapper),
            F32 => dtype_case!(f32, $value, $wrapper),
            F64 => dtype_case!(f64, $value, $wrapper),
            Bool => dtype_case!(bool, $value, $wrapper),
            Char => dtype_case!(char, $value, $wrapper),
            Str => dtype_case!(String, $value, $wrapper),
        }
    }};

    ($dtype:expr, $value:expr, $wrapper:expr; $functor:ident) => {{
        match $dtype {
            USIZE => dtype_case!($functor<usize>, $value, $wrapper),
            U8 => dtype_case!($functor<u8>, $value, $wrapper),
            U16 => dtype_case!($functor<u16>, $value, $wrapper),
            U32 => dtype_case!($functor<u32>, $value, $wrapper),
            U64 => dtype_case!($functor<u64>, $value, $wrapper),
            ISIZE => dtype_case!($functor<isize>, $value, $wrapper),
            I8 => dtype_case!($functor<i8>, $value, $wrapper),
            I16 => dtype_case!($functor<i16>, $value, $wrapper),
            I32 => dtype_case!($functor<i32>, $value, $wrapper),
            I64 => dtype_case!($functor<i64>, $value, $wrapper),
            F32 => dtype_case!($functor<f32>, $value, $wrapper),
            F64 => dtype_case!($functor<f64>, $value, $wrapper),
            Bool => dtype_case!($functor<bool>, $value, $wrapper),
            Char => dtype_case!($functor<char>, $value, $wrapper),
            Str => dtype_case!($functor<String>, $value, $wrapper),
        }
    }};
}

fn len<T>(x: Vec<T>) -> usize {
    x.len()
}

fn to_string<T: fmt::Display>(x: T) -> String {
    x.to_string()
}

// =============================================================================
// Implementations
// =============================================================================

impl Scalar {
    pub fn to_series(self) -> Series {
        dtype_match!(self.dtype, vec![self.unwrap()], Series::new; Vec)
    }

    pub fn to_string(self) -> String {
        dtype_match!(self.dtype, self.unwrap(), to_string)
    }
}

impl Series {
    pub fn at(&self, i: usize) -> Scalar {
        dtype_match!(self.dtype, self.at_raw(i), Scalar::new)
    }

    pub fn len(&self) -> usize {
        dtype_match!(self.dtype, self.as_slice().to_vec(), len; Vec)
    }
}

impl_typed_scalar!(usize, USIZE);
impl_typed_scalar!(u8, U8);
impl_typed_scalar!(u16, U16);
impl_typed_scalar!(u32, U32);
impl_typed_scalar!(u64, U64);
impl_typed_scalar!(isize, ISIZE);
impl_typed_scalar!(i8, I8);
impl_typed_scalar!(i16, I16);
impl_typed_scalar!(i32, I32);
impl_typed_scalar!(i64, I64);
impl_typed_scalar!(f32, F32);
impl_typed_scalar!(f64, F64);
impl_typed_scalar!(bool, Bool);
impl_typed_scalar!(char, Char);
impl_typed_scalar!(String, Str);

impl_typed_vector!(usize, USIZE);
impl_typed_vector!(u8, U8);
impl_typed_vector!(u16, U16);
impl_typed_vector!(u32, U32);
impl_typed_vector!(u64, U64);
impl_typed_vector!(isize, ISIZE);
impl_typed_vector!(i8, I8);
impl_typed_vector!(i16, I16);
impl_typed_vector!(i32, I32);
impl_typed_vector!(i64, I64);
impl_typed_vector!(f32, F32);
impl_typed_vector!(f64, F64);
impl_typed_vector!(bool, Bool);
impl_typed_vector!(char, Char);
impl_typed_vector!(String, Str);

impl DataFrame {
    pub fn new(v: Vec<Series>) -> Self {
        let ics = (0usize .. v.len()).map(|x| x.to_string()).collect();

        Self {
            data: v,
            ics,
        }
    }

    pub fn header(&self) -> &Vec<String> {
        &self.ics
    }

    pub fn header_mut(&mut self) -> &mut Vec<String> {
        &mut self.ics
    }

    pub fn set_header(&mut self, new_header: Vec<&str>) {
        assert_eq!(self.ics.len(), new_header.len(), "Improper Header length!");
        self.ics = new_header.into_iter().map(|x| x.to_string()).collect();
    }

    pub fn push(&mut self, name: &str, series: Series) {
        if self.ics.len() > 0 {
            assert_eq!(self.ics.iter().find(|x| x.as_str() == name), None, "Repetitive index!");
        }
        self.ics.push(name.to_string());
        self.data.push(series);
    }

    pub fn row(&self, i: usize) -> DataFrame {
        let mut df = DataFrame::new(vec![]);
        for (j, series) in self.data.iter().enumerate() {
            let s = series.at(i);
            let new_series = s.to_series();
            df.push(&self.ics[j], new_series);
        }
        df
    }

    pub fn spread(&self) -> String {
        let r: usize = self.data.iter().fold(0, |max_len, column| max(max_len, column.len()));
        let h = self.header();
        
        let mut result = String::new();

        if r > 100 {
            let lc1 = ((r as f64).log10() as usize) + 5;
            result.push_str(&tab("", lc1));

            let mut space_vec: Vec<usize> = vec![];
            for i in 0 .. self.data.len() {
                let v = &self[i];
                let mut space = 0usize;
                for j in 0 .. 5 {
                    let elem = v.at(j);
                    match v.dtype {
                        F32 => {
                            let elem: f32 = elem.unwrap();
                            space = max(
                                space,
                                min(format!("{:.4}", elem).len(), elem.to_string().len())
                            );
                        }
                        F64 => {
                            let elem: f64 = elem.unwrap();
                            space = max(
                                space,
                                min(format!("{:.4}", elem).len(), elem.to_string().len())
                            );
                        }
                        _ => {
                            space = elem.to_string().len();  
                        }
                    } // match
                } // for 0 .. 5
                if v.len() >= r-5 {
                    for j in v.len()-5 .. v.len() {
                        let elem = v.at(j);
                        match v.dtype {
                            F32 => {
                                let elem: f32 = elem.unwrap();
                                space = max(
                                    space,
                                    min(format!("{:.4}", elem).len(), elem.to_string().len())
                                );
                            }
                            F64 => {
                                let elem: f64 = elem.unwrap();
                                space = max(
                                    space,
                                    min(format!("{:.4}", elem).len(), elem.to_string().len())
                                );
                            }
                            _ => {
                                space = elem.to_string().len();  
                            }
                        }
                    }
                } // if v.len >= r-5
                space = max(space + 1, 5);
                let k = &h[i];
                if k.len() >= space {
                    space = k.len() + 1;
                }
                result.push_str(&tab(k, space));
                space_vec.push(space);
            } // for i in 0 .. len()
            result.push('\n');
            
            for i in 0 .. 5 {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for j in 0 .. self.data.len() {
                    let v = &self[j];
                    let space = space_vec[j];
                    if i < v.len() {
                        let elem = v.at(i);
                        let st = match elem.dtype {
                            F32 => {
                                let elem: f32 = elem.unwrap();
                                let st1 = format!("{:.4}", elem);
                                let st2 = elem.to_string();

                                if st1.len() < st2.len() {
                                    st1
                                } else {
                                    st2
                                }
                            }
                            F64 => {
                                let elem: f64 = elem.unwrap();
                                let st1 = format!("{:.4}", elem);
                                let st2 = elem.to_string();

                                if st1.len() < st2.len() {
                                    st1
                                } else {
                                    st2
                                }
                            }
                            _ => {
                                elem.to_string()
                            }
                        };

                        result.push_str(&tab(&st, space));
                    }  else { 
                        result.push_str(&tab("", space));      
                    } // if i < v.len()
                } // for j in 0 .. self.data.len()
                result.push('\n');
            } // for i in 0 .. 5
            result.push_str(&tab("...", lc1));
            for j in 0 .. self.data.len() {
                let space = space_vec[j];
                result.push_str(&tab("...", space));
            } // for j in 0 .. self.data.len()
            result.push('\n');
            for i in r-5 .. r {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for j in 0 .. self.data.len() {
                    let v = &self[j];
                    let space = space_vec[j];
                    if i < v.len() {
                        let elem = v.at(i);
                        let st = match elem.dtype {
                            F32 => {
                                let elem: f32 = elem.unwrap();
                                let st1 = format!("{:.4}", elem);
                                let st2 = elem.to_string();

                                if st1.len() < st2.len() {
                                    st1
                                } else {
                                    st2
                                }
                            }
                            F64 => {
                                let elem: f64 = elem.unwrap();
                                let st1 = format!("{:.4}", elem);
                                let st2 = elem.to_string();

                                if st1.len() < st2.len() {
                                    st1
                                } else {
                                    st2
                                }
                            }
                            _ => {
                                elem.to_string()
                            }
                        };
                        result.push_str(&tab(&st, space));
                    } else {
                        result.push_str(&tab("", space));
                    }
                }
                if i == r-1 {
                    break;
                }
                result.push('\n');
            }
            return result;
        } // if r > 100

        result.push_str(&tab("", 5));
        let mut space_vec: Vec<usize> = vec![];

        for i in 0 .. self.data.len() {
            let v = &self[i];
            let mut space = 0usize;
            for j in 0 .. v.len() {
                let elem = v.at(j);
                match v.dtype {
                    F32 => {
                        let elem: f32 = elem.unwrap();
                        space = max(
                            space,
                            min(format!("{:.4}", elem).len(), elem.to_string().len())
                        );
                    }
                    F64 => {
                        let elem: f64 = elem.unwrap();
                        space = max(
                            space,
                            min(format!("{:.4}", elem).len(), elem.to_string().len())
                        );
                    }
                    _ => {
                        space = elem.to_string().len();  
                    }
                } // match
            }
            space = max(space + 1, 5);
            let k = &h[i];
            if k.len() >= space {
                space = k.len() + 1;
            }
            result.push_str(&tab(k, space));
            space_vec.push(space);
        } // for i in 0 .. len()
        result.push('\n');

        for i in 0 .. r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for j in 0 .. self.data.len() {
                let v = &self[j];
                let space = space_vec[j];
                if i < v.len() {
                    let elem = v.at(i);
                    let st = match elem.dtype {
                        F32 => {
                            let elem: f32 = elem.unwrap();
                            let st1 = format!("{:.4}", elem);
                            let st2 = elem.to_string();

                            if st1.len() < st2.len() {
                                st1
                            } else {
                                st2
                            }
                        }
                        F64 => {
                            let elem: f64 = elem.unwrap();
                            let st1 = format!("{:.4}", elem);
                            let st2 = elem.to_string();

                            if st1.len() < st2.len() {
                                st1
                            } else {
                                st2
                            }
                        }
                        _ => {
                            elem.to_string()
                        }
                    };
                    result.push_str(&tab(&st, space));
                } else {
                    result.push_str(&tab("", space));
                }
            }
            if i == (r - 1) {
                break;
            }
            result.push('\n');
        }
        result
    }
}

impl Index<&str> for DataFrame {
    type Output = Series;

    fn index(&self, index: &str) -> &Self::Output {
        let i = self.ics.iter().position(|x| x.as_str() == index).unwrap();
        &self.data[i]
    }
}

impl IndexMut<&str> for DataFrame {
    fn index_mut(&mut self, index: &str) -> &mut Self::Output {
        let i = self.ics.iter().position(|x| x.as_str() == index).unwrap();
        &mut self.data[i]
    }
}

impl Index<usize> for DataFrame {
    type Output = Series;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl IndexMut<usize> for DataFrame {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl fmt::Display for DataFrame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}
