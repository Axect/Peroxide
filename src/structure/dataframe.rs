//use self::csv::{ReaderBuilder, WriterBuilder};
//use crate::structure::matrix::{matrix, Matrix, Shape::*};
//use crate::traits::fp::FPMatrix;
//use crate::util::useful::tab;
//use indexmap::{map::Keys, IndexMap};
//use json::JsonValue;
//use std::cmp::{max, min};
//use std::collections::HashMap;
//use std::error::Error;
use std::ops::{Index, IndexMut};
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

pub struct IndexSeries {
    pub name: String,
    pub ics: Vec<String>,
}

// =============================================================================
// Implementations
// =============================================================================

impl DTypeValue {
    pub fn to_array(&self) -> DTypeArray {
        match self {
            DTypeValue::USIZE(x) => DTypeArray::USIZE(vec![*x]),
            DTypeValue::ISIZE(x) => DTypeArray::ISIZE(vec![*x]),
            DTypeValue::U8(x) => DTypeArray::U8(vec![*x]),
            DTypeValue::I8(x) => DTypeArray::I8(vec![*x]),
            DTypeValue::U16(x) => DTypeArray::U16(vec![*x]),
            DTypeValue::I16(x) => DTypeArray::I16(vec![*x]),
            DTypeValue::U32(x) => DTypeArray::U32(vec![*x]),
            DTypeValue::I32(x) => DTypeArray::I32(vec![*x]),
            DTypeValue::F32(x) => DTypeArray::F32(vec![*x]),
            DTypeValue::U64(x) => DTypeArray::U64(vec![*x]),
            DTypeValue::I64(x) => DTypeArray::I64(vec![*x]),
            DTypeValue::F64(x) => DTypeArray::F64(vec![*x]),
            DTypeValue::Bool(x) => DTypeArray::Bool(vec![*x]),
            DTypeValue::Char(x) => DTypeArray::Char(vec![*x]),
            DTypeValue::Str(x) => DTypeArray::Str(vec![x.clone()]),
        }
    }

    pub fn unwrap_usize(&self) -> usize {
        match self {
            DTypeValue::USIZE(x) => *x,
            _ => panic!("Can't unwrap usize"),
        }
    }

    pub fn unwrap_u8(&self) -> u8 {
        match self {
            DTypeValue::U8(x) => *x,
            _ => panic!("Can't unwrap u8"),
        }
    }

    pub fn unwrap_u16(&self) -> u16 {
        match self {
            DTypeValue::U16(x) => *x,
            _ => panic!("Can't unwrap u16"),
        }
    }

    pub fn unwrap_u32(&self) -> u32 {
        match self {
            DTypeValue::U32(x) => *x,
            _ => panic!("Can't unwrap u32"),
        }
    }

    pub fn unwrap_u64(&self) -> u64 {
        match self {
            DTypeValue::U64(x) => *x,
            _ => panic!("Can't unwrap u64"),
        }
    }

    pub fn unwrap_isize(&self) -> isize {
        match self {
            DTypeValue::ISIZE(x) => *x,
            _ => panic!("Can't unwrap isize"),
        }
    }

    pub fn unwrap_i8(&self) -> i8 {
        match self {
            DTypeValue::I8(x) => *x,
            _ => panic!("Can't unwrap i8"),
        }
    }

    pub fn unwrap_i16(&self) -> i16 {
        match self {
            DTypeValue::I16(x) => *x,
            _ => panic!("Can't unwrap i16"),
        }
    }

    pub fn unwrap_i32(&self) -> i32 {
        match self {
            DTypeValue::I32(x) => *x,
            _ => panic!("Can't unwrap i32"),
        }
    }

    pub fn unwrap_i64(&self) -> i64 {
        match self {
            DTypeValue::I64(x) => *x,
            _ => panic!("Can't unwrap i64"),
        }
    }

    pub fn unwrap_f32(&self) -> f32 {
        match self {
            DTypeValue::F32(x) => *x,
            _ => panic!("Can't unwrap f32"),
        }
    }

    pub fn unwrap_f64(&self) -> f64 {
        match self {
            DTypeValue::F64(x) => *x,
            _ => panic!("Can't unwrap f64"),
        }
    }

    pub fn unwrap_bool(&self) -> bool {
        match self {
            DTypeValue::Bool(x) => *x,
            _ => panic!("Can't unwrap bool"),
        }
    }

    pub fn unwrap_str(&self) -> String {
        match self {
            DTypeValue::Str(x) => x.clone(),
            _ => panic!("Can't unwrap String"),
        }
    }

    pub fn unwrap_char(&self) -> char {
        match self {
            DTypeValue::Char(x) => *x,
            _ => panic!("Can't unwrap char"),
        }
    }
}

impl Series {
    pub fn new(values: DTypeArray) -> Self {
        let dtype = match values {
            DTypeArray::USIZE(_) => DType::USIZE,
            DTypeArray::U8(_) => DType::U8,
            DTypeArray::U16(_) => DType::U16,
            DTypeArray::U32(_) => DType::U32,
            DTypeArray::U64(_) => DType::U64,
            DTypeArray::ISIZE(_) => DType::ISIZE,
            DTypeArray::I8(_) => DType::I8,
            DTypeArray::I16(_) => DType::I16,
            DTypeArray::I32(_) => DType::I32,
            DTypeArray::I64(_) => DType::I64,
            DTypeArray::F32(_) => DType::F32,
            DTypeArray::F64(_) => DType::F64,
            DTypeArray::Bool(_) => DType::Bool,
            DTypeArray::Str(_) => DType::Str,
            DTypeArray::Char(_) => DType::Char,
        };

        Self {
            values,
            dtype,
        }
    }

    pub fn at(&self, i: usize) -> DTypeValue {
        match &self.values {
            DTypeArray::USIZE(v) => DTypeValue::USIZE(v[i]),
            DTypeArray::U8(v) => DTypeValue::U8(v[i]),
            DTypeArray::U16(v) => DTypeValue::U16(v[i]),
            DTypeArray::U32(v) => DTypeValue::U32(v[i]),
            DTypeArray::U64(v) => DTypeValue::U64(v[i]),
            DTypeArray::ISIZE(v) => DTypeValue::ISIZE(v[i]),
            DTypeArray::I8(v) => DTypeValue::I8(v[i]),
            DTypeArray::I16(v) => DTypeValue::I16(v[i]),
            DTypeArray::I32(v) => DTypeValue::I32(v[i]),
            DTypeArray::I64(v) => DTypeValue::I64(v[i]),
            DTypeArray::F32(v) => DTypeValue::F32(v[i]),
            DTypeArray::F64(v) => DTypeValue::F64(v[i]),
            DTypeArray::Bool(v) => DTypeValue::Bool(v[i]),
            DTypeArray::Str(v) => DTypeValue::Str(v[i].clone()),
            DTypeArray::Char(v) => DTypeValue::Char(v[i]),
        }
    }
}

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
            let new_series = Series::new(series.at(i).to_array());
            df.push(&self.ics[j], new_series);
        }
        df
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
