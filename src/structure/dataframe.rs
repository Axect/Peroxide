//use self::csv::{ReaderBuilder, WriterBuilder};
//use crate::structure::matrix::{matrix, Matrix, Shape::*};
//use crate::traits::fp::FPMatrix;
//use indexmap::{map::Keys, IndexMap};
//use json::JsonValue;
use std::collections::HashMap;
use std::fmt;
use std::ops::{Index, IndexMut};
use std::cmp::{max, min};
use std::error::Error;
use crate::util::useful::tab;
use csv::{ReaderBuilder, WriterBuilder};
use DType::{
    USIZE,U8,U16,U32,U64,
    ISIZE,I8,I16,I32,I64,
    F32,F64,Bool,Char,Str
};

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

// =============================================================================
// Traits
// =============================================================================
pub trait TypedScalar<T> {
    fn new(s: T) -> Self where Self: Sized;
    fn unwrap(self) -> T;
}

pub trait TypedVector<T> {
    fn new(v: Vec<T>) -> Self;
    fn to_vec(&self) -> Vec<T>;
    fn as_slice(&self) -> &[T];
    fn as_slice_mut(&mut self) -> &mut [T];
    fn at_raw(&self, i: usize) -> T;
    fn push(&mut self, elem: T);
}

// =============================================================================
// Macros & Private functions
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

            fn to_vec(&self) -> Vec<$type> {
                self.as_slice().to_vec()
            }

            fn as_slice(&self) -> &[$type] {
                match &self.values {
                    DTypeArray::$dtype(v) => v,
                    _ => panic!("Can't convert to {:?} vector", $dtype),
                }
            }

            fn as_slice_mut(&mut self) -> &mut [$type] {
                match &mut self.values {
                    DTypeArray::$dtype(v) => v,
                    _ => panic!("Can't convert to {:?} vector", $dtype),
                }
            }

            fn at_raw(&self, i: usize) -> $type {
                let v: &[$type] = self.as_slice();
                v[i].clone()
            }

            fn push(&mut self, elem: $type) {
                let v: &mut Vec<$type> = match &mut self.values {
                    DTypeArray::$dtype(v) => v,
                    _ => panic!("Can't convert to {:?} vector", $dtype),
                };
                v.push(elem);
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

macro_rules! set_space {
    ($elem:expr) => {{
        match $elem.dtype {
            F32 => {
                let elem: f32 = $elem.unwrap();
                let st1 = format!("{:.4}", elem);
                let st2 = elem.to_string();

                if st1.len() < st2.len() {
                    st1
                } else {
                    st2
                }
            }
            F64 => {
                let elem: f64 = $elem.unwrap();
                let st1 = format!("{:.4}", elem);
                let st2 = elem.to_string();

                if st1.len() < st2.len() {
                    st1
                } else {
                    st2
                }
            }
            _ => $elem.to_string()
        }
    }};

    ($elem:expr, $space:expr) => {{
        match $elem.dtype {
            F32 => {
                let elem: f32 = $elem.unwrap();
                $space = max(
                    $space,
                    min(format!("{:.4}", elem).len(), elem.to_string().len())
                );
            }
            F64 => {
                let elem: f64 = $elem.unwrap();
                $space = max(
                    $space,
                    min(format!("{:.4}", elem).len(), elem.to_string().len())
                );
            }
            _ => {
                $space = $elem.to_string().len();
            }
        }
    }};
}

macro_rules! format_float_vec {
    ($self:expr) => {{
        let mut result = String::new();
        result.push_str("[");
        for i in 0 .. $self.len() {
            let st1 = format!("{:.4}", $self[i]);
            let st2 = $self[i].to_string();
            let st = if st1.len() < st2.len() {
                st1
            } else {
                st2
            };
            result.push_str(&st);
            if i == $self.len() - 1 {
                break;
            }
            result.push_str(", ");
        }
        result.push_str("]");
        result
    }};
}

/// ty1 -> ty2
macro_rules! type_cast_vec {
    ($ty1:ty, $ty2:ty, $to_vec:expr, $wrapper:expr) => {{
        let y: Vec<$ty1> = $to_vec;
        let x: Vec<$ty2> = y.into_iter().map(|x| x as $ty2).collect();
        $wrapper(x)
    }};
}

macro_rules! string_cast_vec {
    ($ty1:ty, $to_vec:expr, $wrapper:expr) => {{
        let y: Vec<$ty1> = $to_vec;
        let x: Vec<String> = y.into_iter().map(|x| x.to_string()).collect();
        $wrapper(x)
    }}
}

macro_rules! type_parse_vec {
    ($ty2:ty, $to_vec:expr, $wrapper:expr) => {{
        let y: Vec<String> = $to_vec.to_vec();
        let x: Vec<$ty2> = y.into_iter().map(|x| x.parse().unwrap()).collect();
        $wrapper(x)
    }};
}

macro_rules! dtype_parse_vec_part {
    ($dt2:expr, $to_vec:expr, $wrapper:expr) => {{
        match $dt2 {
            USIZE => type_parse_vec!(usize, $to_vec, $wrapper),
            U8 => type_parse_vec!(u8, $to_vec, $wrapper),
            U16 => type_parse_vec!(u16, $to_vec, $wrapper),
            U32 => type_parse_vec!(u32, $to_vec, $wrapper),
            U64 => type_parse_vec!(u64, $to_vec, $wrapper),
            ISIZE => type_parse_vec!(isize, $to_vec, $wrapper),
            I8 => type_parse_vec!(i8, $to_vec, $wrapper),
            I16 => type_parse_vec!(i16, $to_vec, $wrapper),
            I32 => type_parse_vec!(i32, $to_vec, $wrapper),
            I64 => type_parse_vec!(i64, $to_vec, $wrapper),
            F32 => type_parse_vec!(f32, $to_vec, $wrapper),
            F64 => type_parse_vec!(f64, $to_vec, $wrapper),
            Bool => type_parse_vec!(bool, $to_vec, $wrapper),
            Char => type_parse_vec!(char, $to_vec, $wrapper),
            Str => type_parse_vec!(String, $to_vec, $wrapper),
        }
    }};
}

macro_rules! dtype_cast_vec_part {
    ($ty1:ty, $dt2:expr, $to_vec:expr, $wrapper:expr) => {{
        match $dt2 {
            USIZE => type_cast_vec!($ty1, usize, $to_vec, $wrapper),
            U8 => type_cast_vec!($ty1, u8, $to_vec, $wrapper),
            U16 => type_cast_vec!($ty1, u16, $to_vec, $wrapper),
            U32 => type_cast_vec!($ty1, u32, $to_vec, $wrapper),
            U64 => type_cast_vec!($ty1, u64, $to_vec, $wrapper),
            ISIZE => type_cast_vec!($ty1, isize, $to_vec, $wrapper),
            I8 => type_cast_vec!($ty1, i8, $to_vec, $wrapper),
            I16 => type_cast_vec!($ty1, i16, $to_vec, $wrapper),
            I32 => type_cast_vec!($ty1, i32, $to_vec, $wrapper),
            I64 => type_cast_vec!($ty1, i64, $to_vec, $wrapper),
            F32 => type_cast_vec!($ty1, f32, $to_vec, $wrapper),
            F64 => type_cast_vec!($ty1, f64, $to_vec, $wrapper),
            Str => string_cast_vec!($ty1, $to_vec, $wrapper),
            _ => panic!("Can't convert type to {}", $dt2),
        }
    }};
}

macro_rules! dtype_cast_vec {
    ($dt1:expr, $dt2:expr, $to_vec:expr, $wrapper:expr) => {{
        match $dt1 {
            USIZE => dtype_cast_vec_part!(usize, $dt2, $to_vec, $wrapper),
            U8 => dtype_cast_vec_part!(u8, $dt2, $to_vec, $wrapper),
            U16 => dtype_cast_vec_part!(u16, $dt2, $to_vec, $wrapper),
            U32 => dtype_cast_vec_part!(u32, $dt2, $to_vec, $wrapper),
            U64 => dtype_cast_vec_part!(u64, $dt2, $to_vec, $wrapper),
            ISIZE => dtype_cast_vec_part!(isize, $dt2, $to_vec, $wrapper),
            I8 => dtype_cast_vec_part!(i8, $dt2, $to_vec, $wrapper),
            I16 => dtype_cast_vec_part!(i16, $dt2, $to_vec, $wrapper),
            I32 => dtype_cast_vec_part!(i32, $dt2, $to_vec, $wrapper),
            I64 => dtype_cast_vec_part!(i64, $dt2, $to_vec, $wrapper),
            F32 => dtype_cast_vec_part!(f32, $dt2, $to_vec, $wrapper),
            F64 => dtype_cast_vec_part!(f64, $dt2, $to_vec, $wrapper),
            Str => dtype_parse_vec_part!($dt2, $to_vec, $wrapper),
            Char => {
                match $dt2 {
                    Str => string_cast_vec!(char, $to_vec, $wrapper),
                    _ => panic!("Can't convert char type to {}", $dt2),
                }
            }
            Bool => {
                match $dt2 {
                    Str => string_cast_vec!(bool, $to_vec, $wrapper),
                    _ => panic!("Can't convert bool type to {}", $dt2),
                }
            }
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
// Implementations of DType variables
// =============================================================================
impl DType {
    pub fn is_numeric(&self) -> bool {
        match self {
            Bool => false,
            Str => false,
            Char => false,
            _ => true,
        }
    }

    pub fn is_integer(&self) -> bool {
        match self {
            Bool => false,
            Str => false,
            Char => false,
            F32 => false,
            F64 => false,
            _ => true,
        }
    }
}

impl fmt::Display for DType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let st = match self {
            USIZE => "usize",
            U8 => "u8",
            U16 => "u16",
            U32 => "u32",
            U64 => "u64",
            ISIZE => "isize",
            I8 => "i8",
            I16 => "i16",
            I32 => "i32",
            I64 => "i64",
            F32 => "f32",
            F64 => "f64",
            Bool => "bool",
            Char => "char",
            Str => "String",
        };
        write!(f, "{}", st)
    }
}

impl fmt::Display for DTypeArray {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let st = match self {
            DTypeArray::USIZE(v) => format!("array: {:?}\ndtype: usize", v),
            DTypeArray::U8(v) => format!("array: {:?}\ndtype: u8", v),
            DTypeArray::U16(v) => format!("array: {:?}\ndtype: u16", v),
            DTypeArray::U32(v) => format!("array: {:?}\ndtype: u32", v),
            DTypeArray::U64(v) => format!("array: {:?}\ndtype: u64", v),
            DTypeArray::ISIZE(v) => format!("array: {:?}\ndtype: isize", v),
            DTypeArray::I8(v) => format!("array: {:?}\ndtype: i8", v),
            DTypeArray::I16(v) => format!("array: {:?}\ndtype: i16", v),
            DTypeArray::I32(v) => format!("array: {:?}\ndtype: i32", v),
            DTypeArray::I64(v) => format!("array: {:?}\ndtype: i64", v),
            DTypeArray::F32(v) => format!("array: {}\ndtype: f32", format_float_vec!(v)),
            DTypeArray::F64(v) => format!("array: {}\ndtype: f64", format_float_vec!(v)),
            DTypeArray::Bool(v) => format!("array: {:?}\ndtype: bool", v),
            DTypeArray::Str(v) => format!("array: {:?}\ndtype: String", v),
            DTypeArray::Char(v) => format!("array: {:?}\ndtype: char", v),
        };
        write!(f, "{}", st)
    }
}

// =============================================================================
// Implementations for Scalar & Series
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

    pub fn to_type(&self, dtype: DType) -> Series {
        dtype_cast_vec!(self.dtype, dtype, self.to_vec(), Series::new)
    }

    pub fn as_type(&mut self, dtype: DType) {
        let x = self.to_type(dtype);
        self.dtype = x.dtype;
        self.values = x.values;
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

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let st = format!("{}, dtype:{}", self.clone().to_string(), self.dtype);
        write!(f, "{}", st)
    }
}

// =============================================================================
// Implementation for DataFrame
// =============================================================================

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
                    set_space!(elem, space);
                }
                if v.len() >= r-5 {
                    for j in v.len()-5 .. v.len() {
                        let elem = v.at(j);
                        set_space!(elem, space);
                    }
                }
                space = max(space + 1, 5);
                let k = &h[i];
                if k.len() >= space {
                    space = k.len() + 1;
                }
                result.push_str(&tab(k, space));
                space_vec.push(space);
            }
            result.push('\n');
            
            for i in 0 .. 5 {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for j in 0 .. self.data.len() {
                    let v = &self[j];
                    let space = space_vec[j];
                    if i < v.len() {
                        let elem = v.at(i);
                        let st = set_space!(elem);
                        result.push_str(&tab(&st, space));
                    }  else { 
                        result.push_str(&tab("", space));      
                    }
                }
                result.push('\n');
            }
            result.push_str(&tab("...", lc1));
            for j in 0 .. self.data.len() {
                let space = space_vec[j];
                result.push_str(&tab("...", space));
            }
            result.push('\n');
            for i in r-5 .. r {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for j in 0 .. self.data.len() {
                    let v = &self[j];
                    let space = space_vec[j];
                    if i < v.len() {
                        let elem = v.at(i);
                        let st = set_space!(elem);
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
        }

        result.push_str(&tab("", 5));
        let mut space_vec: Vec<usize> = vec![];

        for i in 0 .. self.data.len() {
            let v = &self[i];
            let mut space = 0usize;
            for j in 0 .. v.len() {
                let elem = v.at(j);
                set_space!(elem, space)
            }
            space = max(space + 1, 5);
            let k = &h[i];
            if k.len() >= space {
                space = k.len() + 1;
            }
            result.push_str(&tab(k, space));
            space_vec.push(space);
        }
        result.push('\n');

        for i in 0 .. r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for j in 0 .. self.data.len() {
                let v = &self[j];
                let space = space_vec[j];
                if i < v.len() {
                    let elem = v.at(i);
                    let st = set_space!(elem);
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

    pub fn as_types(&mut self, dtypes: Vec<DType>) {
        assert_eq!(self.data.len(), dtypes.len(), "Length of dtypes are not compatible with DataFrame");
        for (i, dtype) in dtypes.into_iter().enumerate() {
            self[i].as_type(dtype);
        }
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

// =============================================================================
// IO Implementations
// =============================================================================

pub trait WithCSV: Sized {
    fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>>;
    fn read_csv(file_path: &str, delimiter: char) -> Result<Self, Box<dyn Error>>;
}

impl WithCSV for DataFrame {
    fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r: usize = self
            .data
            .iter()
            .fold(0, |max_len, column| max(max_len, column.len()));
        let c: usize = self.data.len();
        wtr.write_record(
            self.header().clone()
        )?;
        
        for i in 0 .. r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for (j, v) in self.data.iter().enumerate() {
                if i < v.len() {
                    record[j] = v.at(i).to_string();
                }
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    fn read_csv(file_path: &str, delimiter: char) -> Result<Self, Box<dyn Error>> {
        let mut rdr = ReaderBuilder::new()
            .has_headers(true)
            .delimiter(delimiter as u8)
            .from_path(file_path)?;

        let headers_vec = rdr.headers()?;
        let headers = headers_vec.iter().map(|x| x).collect::<Vec<&str>>();
        let mut result = DataFrame::new(vec![]);
        for h in headers.iter() {
            result.push(*h, Series::new(Vec::<String>::new()));
        }

        for rec in rdr.deserialize() {
            let record: HashMap<String, String> = rec?;
            for head in record.keys() {
                let value = &record[head];
                if value.len() > 0 {
                    result[head.as_str()].push(value.to_string());
                }
            }
        }

        Ok(result)
    }
}
