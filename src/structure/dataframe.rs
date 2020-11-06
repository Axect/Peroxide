extern crate csv;

use self::csv::{ReaderBuilder, WriterBuilder};
use crate::structure::matrix::{matrix, Matrix, Shape::*};
use crate::traits::fp::FPMatrix;
use crate::util::useful::tab;
use indexmap::{map::Keys, IndexMap};
use json::JsonValue;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::error::Error;
use std::ops::{Index, IndexMut};
use std::{fmt, fmt::Debug};

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

pub struct Series {
    pub name: String,
    pub index: IndexSeries,
    pub values: DTypeArray
}

pub struct IndexSeries {
    pub name: String,
    pub ics: Vec<String>,
}