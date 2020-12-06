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

    pub fn header_replace(&mut self, new_header: Vec<&str>) {
        assert_eq!(self.ics.len(), new_header.len(), "Improper Header length!");
        self.ics = new_header.into_iter().map(|x| x.to_string()).collect();
    }

    pub fn push(&mut self, name: &str, series: Series) {
        assert_eq!(self.ics.iter().find(|x| x.as_str() == name), None, "Repetitive index!");
        
        self.ics.push(name.to_string());
        self.data.push(series);
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