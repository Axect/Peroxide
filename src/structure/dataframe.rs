//! Pandas-like dataframe & series.
//!
//! ## Series
//!
//! ### 1. Declare Series
//!
//! * To declare series, you should have `Vec<T>` where `T` is one of following types.
//!
//! | Primitive type | DType   |
//! | :-----:  | :-----: |
//! | `usize`  | `USIZE` |
//! | `u8`     | `U8`    |
//! | `u16`    | `U16`   |
//! | `u32`    | `U32`   |
//! | `u64`    | `U64`   |
//! | `isize`  | `ISIZE` |
//! | `i8`     | `I8`    |
//! | `i16`    | `I16`   |
//! | `i32`    | `I32`   |
//! | `i64`    | `I64`   |
//! | `f32`    | `F32`   |
//! | `f64`    | `F64`   |
//! | `bool`   | `Bool`  |
//! | `char`   | `Char`  |
//! | `String` | `Str`   |
//!
//! * If you prepare `Vec<T>`, then `Series::new(Vec<T>)`
//!
//! ### 2. Methods for Series
//!
//! * `TypedVector<T> trait for Series`
//!     
//!     ```ignore
//!     pub trait TypedVector<T> {
//!         fn new(v: Vec<T>) -> Self;
//!         fn to_vec(&self) -> Vec<T>;
//!         fn as_slice(&self) -> &[T];
//!         fn as_slice_mut(&mut self) -> &mut [T];
//!         fn at_raw(&self, i: usize) -> T;
//!         fn push(&mut self, elem: T);
//!     }
//!     ```
//!
//! * `Series` methods
//!
//!     ```ignore
//!     impl Series {
//!         // Core
//!         pub fn at(&self, i: usize) -> Scalar;
//!         pub fn len(&self) -> usize;
//!         pub fn to_type(&self, dtype: DType) -> Series;
//!         pub fn as_type(&mut self, dtype: DType);
//!         pub fn select_indices(&self, indices: &[usize]) -> Series;
//!         pub fn to_f64_vec(&self) -> anyhow::Result<Vec<f64>>;
//!
//!         // Statistics (numeric types only, except min/max)
//!         pub fn sum(&self) -> anyhow::Result<f64>;
//!         pub fn mean(&self) -> anyhow::Result<f64>;
//!         pub fn var(&self) -> anyhow::Result<f64>;
//!         pub fn sd(&self) -> anyhow::Result<f64>;
//!         pub fn min(&self) -> anyhow::Result<Scalar>;
//!         pub fn max(&self) -> anyhow::Result<Scalar>;
//!     }
//!     ```
//!
//!     * `at` is simple getter for `Series`. It returns `Scalar`.
//!     * `as_type` is a method for mutable type casting.
//!         * All types can be changed to `Str`.
//!         * All integer & float types can be exchanged.
//!         * `Bool, Char` can be changed to `Str` or `U8` only.
//!         * `U8` can be changed to all types.
//!     * `select_indices` selects elements by indices, returning a new Series.
//!     * `to_f64_vec` converts numeric Series to `Vec<f64>` (bridge for statistics).
//!     * `sum`, `mean`, `var`, `sd` convert to `f64` internally via `to_f64_vec`.
//!     * `min`, `max` preserve the original type and return `Scalar`. Works on all ordered types including `Char` and `String`.
//!
//! ### 3. Example
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let a = Series::new(vec![1, 2, 3, 4]);
//!     let b = Series::new(vec!['a', 'b', 'c', 'd']);
//!     let mut c = Series::new(vec![true, false, false, true]);
//!
//!     a.print();       // print for Series
//!     b.dtype.print(); // print for dtype of Series (=Char)
//!     c.as_type(U8);   // Bool => U8
//!
//!     assert_eq!(c.dtype, U8);
//!
//!     // Select by indices
//!     let d = a.select_indices(&[0, 2]);
//!     assert_eq!(d, Series::new(vec![1, 3]));
//!
//!     // Statistics
//!     let e = Series::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
//!     assert_eq!(e.sum().unwrap(), 15.0);
//!     assert_eq!(e.mean().unwrap(), 3.0);
//!     assert_eq!(e.min().unwrap(), Scalar::new(1.0f64));
//!     assert_eq!(e.max().unwrap(), Scalar::new(5.0f64));
//! }
//! ```
//!
//! ## DataFrame
//!
//! ### 1. Declare DataFrame
//!
//! * To declare dataframe, use constructor.
//!     * `DataFrame::new(Vec<Series>)`
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // 1-1. Empty DataFrame
//!     let mut df = DataFrame::new(vec![]);
//!
//!     // 1-2. Push Series
//!     df.push("a", Series::new(vec![1, 2, 3, 4]));
//!     df.push("b", Series::new(vec![0.1, 0.2, 0.3, 0.4]));
//!     df.push("c", Series::new(vec!['a', 'b', 'c', 'd']));
//!
//!     // 1-3. Print
//!     df.print();
//!
//!     // 2-1. Construct Series first
//!     let a = Series::new(vec![1, 2, 3, 4]);
//!     let b = Series::new(vec![0.1, 0.2, 0.3, 0.4]);
//!     let c = Series::new(vec!['a', 'b', 'c', 'd']);
//!
//!     // 2-2. Declare DataFrame with exist Series
//!     let mut dg = DataFrame::new(vec![a, b, c]);
//!
//!     // 2-3. Print or Set header
//!     dg.print();                         // But header: 0 1 2
//!     dg.set_header(vec!["a", "b", "c"]); // Change header
//! }
//! ```
//!
//! ### 2. Methods for DataFrame
//!
//! * `DataFrame` method
//!
//!     ```ignore
//!     impl DataFrame {
//!         // Constructor & Basic
//!         pub fn new(v: Vec<Series>) -> Self;
//!         pub fn header(&self) -> &Vec<String>;
//!         pub fn header_mut(&mut self) -> &mut Vec<String>;
//!         pub fn set_header(&mut self, new_header: Vec<&str>);
//!         pub fn push(&mut self, name: &str, series: Series);
//!         pub fn drop(&mut self, col_header: &str);
//!         pub fn row(&self, i: usize) -> DataFrame;
//!         pub fn spread(&self) -> String;
//!         pub fn as_types(&mut self, dtypes: Vec<DType>);
//!         pub fn filter_by<F>(&self, column: &str, f: F) -> anyhow::Result<DataFrame>;
//!         pub fn mask(&self, mask: &Series) -> anyhow::Result<DataFrame>;
//!         pub fn select_rows(&self, indices: &[usize]) -> DataFrame;
//!
//!         // Shape & Info
//!         pub fn nrow(&self) -> usize;
//!         pub fn ncol(&self) -> usize;
//!         pub fn shape(&self) -> (usize, usize);
//!         pub fn dtypes(&self) -> Vec<DType>;
//!         pub fn is_empty(&self) -> bool;
//!         pub fn contains(&self, col_header: &str) -> bool;
//!
//!         // Row Operations
//!         pub fn head(&self, n: usize) -> DataFrame;
//!         pub fn tail(&self, n: usize) -> DataFrame;
//!         pub fn slice(&self, offset: usize, length: usize) -> DataFrame;
//!
//!         // Column Operations
//!         pub fn select(&self, columns: &[&str]) -> DataFrame;
//!         pub fn rename(&mut self, old: &str, new: &str);
//!         pub fn column_names(&self) -> Vec<&str>;
//!         pub fn select_dtypes(&self, dtypes: &[DType]) -> DataFrame;
//!
//!         // Statistics (numeric columns only)
//!         pub fn describe(&self) -> DataFrame;
//!         pub fn sum(&self) -> DataFrame;
//!         pub fn mean(&self) -> DataFrame;
//!     }
//!     ```
//!
//!     * `push(&mut self, name: &str, series: Series)`: push head & Series pair
//!     * `drop(&mut self, col_header: &str)`: drop specific column by header
//!     * `row(&self, i: usize) -> DataFrame` : Extract $i$-th row as new DataFrame
//!     * `filter_by(&self, column, f)` : Filter DataFrame by specific column
//!     * `mask(&self, mask: &Series)` : Mask DataFrame by boolean Series
//!     * `select_rows(&self, indices)` : Select rows by indices
//!     * `nrow`, `ncol`, `shape` : Row count (max column length), column count, `(nrow, ncol)` tuple
//!     * `dtypes` : `Vec<DType>` of each column's type
//!     * `is_empty` : `true` if no columns or no rows
//!     * `contains(col_header)` : `true` if the column exists
//!     * `head(n)`, `tail(n)` : First / last `n` rows
//!     * `slice(offset, length)` : Row slice starting at `offset`
//!     * `select(columns)` : Select columns by name (panics on missing)
//!     * `rename(old, new)` : Rename a column in-place
//!     * `column_names` : `Vec<&str>` of all headers
//!     * `select_dtypes(dtypes)` : Select columns matching given DTypes
//!     * `describe` : Computes count / mean / sd / min / max for each numeric column
//!     * `sum`, `mean` : Single-row DataFrame with column-wise sum / mean
//!
//! * `WithCSV` trait
//!
//!     ```ignore
//!     pub trait WithCSV: Sized {
//!         fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>>;
//!         fn read_csv(file_path: &str, delimiter: char) -> Result<Self, Box<dyn Error>>;
//!     }
//!     ```
//!
//!     * `csv` feature should be required
//!
//!     ```rust
//!     // Example for CSV
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() -> Result<(), Box<dyn Error>> {
//!     # #[cfg(feature="csv")]
//!     # {
//!         // Write CSV
//!         let mut df = DataFrame::new(vec![]);
//!         df.push("a", Series::new(vec!['x', 'y', 'z']));
//!         df.push("b", Series::new(vec![0, 1, 2]));
//!         df.push("c", Series::new(c!(0.1, 0.2, 0.3)));
//!         df.write_csv("example_data/doc_csv.csv")?;
//!
//!         // Read CSV
//!         let mut dg = DataFrame::read_csv("example_data/doc_csv.csv", ',')?;
//!         dg.as_types(vec![Char, I32, F64]);
//!
//!         assert_eq!(df, dg);
//!     # }
//!
//!         Ok(())
//!     }
//!     ```
//!
//! * `WithNetCDF` trait
//!
//!     ```ignore
//!     pub trait WithNetCDF: Sized {
//!         fn write_nc(&self, file_path: &str) -> Result<(), Box<dyn Error>>;
//!         fn read_nc(file_path: &str) -> Result<Self, Box<dyn Error>>;
//!         fn read_nc_by_header(file_path: &str, header: Vec<&str>) -> Result<Self, Box<dyn Error>>;
//!     }
//!     ```
//!
//!     * `nc` feature should be required
//!     * `libnetcdf` dependency should be required
//!     * `Char`, `Bool` are saved as `U8` type. Thus, for reading `Char` or `Bool` type nc file, explicit type casting is required.
//!
//!     ```
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!
//!     fn main() -> Result<(), Box<dyn Error>> {
//!     #    #[cfg(feature = "nc")]
//!     #    {
//!         // Write netcdf
//!         let mut df = DataFrame::new(vec![]);
//!         df.push("a", Series::new(vec!['x', 'y', 'z']));
//!         df.push("b", Series::new(vec![0, 1, 2]));
//!         df.push("c", Series::new(c!(0.1, 0.2, 0.3)));
//!         df.write_nc("example_data/doc_nc.nc")?;
//!
//!         // Read netcdf
//!         let mut dg = DataFrame::read_nc("example_data/doc_nc.nc")?;
//!         dg["a"].as_type(Char); // Char, Bool are only read/written as U8 type
//!
//!         assert_eq!(df, dg);
//!     #    }
//!
//!         Ok(())
//!     }
//!     ```
//!
//! * `WithParquet` trait
//!
//!     ```ignore
//!     pub trait WithParquet: Sized {
//!         fn write_parquet(&self, file_path: &str, compression: Compression) -> Result<(), Box<dyn Error>>;
//!         fn read_parquet(file_path: &str) -> Result<Self, Box<dyn Error>>;
//!     }
//!     ```
//!
//!     * `parquet` feature should be required
//!     * `Char` is saved with `String` type. Thus, for reading `Char` type parquet file, the output type is `String`.
//!     * **Caution** : For different length `Bool` type column, missing values are filled with `false`.
//!     ```
//!     #[macro_use]
//!     extern crate peroxide;
//!     use peroxide::fuga::*;
//!     
//!     fn main() -> Result<(), Box<dyn Error>> {
//!     #    #[cfg(feature = "parquet")]
//!     #    {
//!         // Write parquet
//!         let mut df = DataFrame::new(vec![]);
//!         df.push("a", Series::new(vec!['x', 'y', 'z']));
//!         df.push("b", Series::new(vec![0, 1, 2]));
//!         df.push("c", Series::new(c!(0.1, 0.2, 0.3)));
//!         df.write_parquet("example_data/doc_pq.parquet", SNAPPY)?;
//!
//!         // Read parquet
//!         let mut dg = DataFrame::read_parquet("example_data/doc_pq.parquet")?;
//!         dg["a"].as_type(Char); // Char is only read/written as String type
//!
//!         assert_eq!(df, dg);
//!     #    }
//!
//!         Ok(())
//!     }
//!     ```

use crate::traits::math::Vector;
use crate::util::{print::LowerExpWithPlus, useful::tab};
#[cfg(feature = "parquet")]
use arrow::datatypes::{
    Float32Type, Float64Type, Int16Type, Int32Type, Int64Type, Int8Type, UInt16Type, UInt32Type,
    UInt64Type, UInt8Type,
};
use std::cmp::{max, min};
#[cfg(feature = "csv")]
use std::collections::HashMap;
#[cfg(feature = "parquet")]
use indexmap::IndexMap;
#[cfg(any(feature = "csv", feature = "nc", feature = "parquet"))]
use std::error::Error;
use std::fmt;
use std::ops::{Index, IndexMut};
#[cfg(feature = "parquet")]
use std::sync::Arc;
use DType::{Bool, Char, Str, F32, F64, I16, I32, I64, I8, ISIZE, U16, U32, U64, U8, USIZE};

#[cfg(feature = "parquet")]
use arrow::{
    array::{Array, BooleanArray, PrimitiveArray, StringArray},
    datatypes::{DataType, Field, Schema},
};
#[cfg(feature = "csv")]
use csv::{ReaderBuilder, WriterBuilder};
#[cfg(feature = "nc")]
use netcdf::{
    types::VariableType,
    variable::{Variable, VariableMut},
    Numeric,
};
#[cfg(feature = "parquet")]
use parquet::{
    arrow::arrow_reader::ParquetRecordBatchReaderBuilder,
    arrow::arrow_writer::compute_leaves,
    arrow::arrow_writer::get_column_writers,
    arrow::arrow_writer::ArrowLeafColumn,
    arrow::ArrowSchemaConverter,
    basic::Compression,
    file::properties::WriterProperties,
    file::writer::{SerializedFileWriter, SerializedRowGroupWriter},
};

// =============================================================================
// Enums
// =============================================================================

/// Data Type enum
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

/// Vector with `DType`
#[derive(Debug, Clone, PartialEq)]
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

/// Scalar with `DType`
#[derive(Debug, Clone, PartialEq, PartialOrd)]
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

/// Generic `DataFrame` structure
///
/// # Example
///
/// ```rust
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     // 1. Series to DataFrame
///     // 1-1. Declare Series
///     let a = Series::new(vec![1, 2, 3, 4]);
///     let b = Series::new(vec![true, false, false, true]);
///     let c = Series::new(vec![0.1, 0.2, 0.3, 0.4]);
///
///     // 1-2. Declare DataFrame (default header: 0, 1, 2)
///     let mut df = DataFrame::new(vec![a, b, c]);
///     df.set_header(vec!["a", "b", "c"]);
///     df.print(); // Pretty print for DataFrame
///
///     // 2. Empty DataFrame
///     let mut dg = DataFrame::new(vec![]);
///     dg.push("a", Series::new(vec![1,2,3,4]));
///     dg.push("b", Series::new(vec![true, false, false, true]));
///     dg.push("c", Series::new(vec![0.1, 0.2, 0.3, 0.4]));
///     dg.print();
///
///     assert_eq!(df, dg);
/// }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct DataFrame {
    pub data: Vec<Series>,
    pub ics: Vec<String>,
}

/// Generic Series
///
/// # Example
///
/// ```rust
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     // Declare Series with Vec<T> (T: primitive type)
///     let a = Series::new(vec![1i32, 2, 3, 4]);
///     a.print();                      // print for Series
///     a.dtype.print();              // print for dtype of Series
///
///     let b: &[i32] = a.as_slice();   // Borrow series to &[T]
///     let c: Vec<i32> = a.to_vec();   // Series to Vec<T> (clone)
///     
///     // ...
/// }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Series {
    pub values: DTypeArray,
    pub dtype: DType,
}

/// Generic Scalar
#[derive(Debug, Clone, PartialEq)]
pub struct Scalar {
    pub value: DTypeValue,
    pub dtype: DType,
}

// =============================================================================
// Traits
// =============================================================================
pub trait TypedScalar<T> {
    fn new(s: T) -> Self
    where
        Self: Sized;
    fn unwrap(self) -> T;
}

pub trait TypedVector<T> {
    fn new(v: Vec<T>) -> Self;
    fn to_vec(&self) -> Vec<T>;
    fn as_slice(&self) -> &[T];
    fn as_slice_mut(&mut self) -> &mut [T];
    fn at_raw(&self, i: usize) -> T;
    fn push(&mut self, elem: T);
    fn map<F: Fn(T) -> T>(&self, f: F) -> Self;
    fn mut_map<F: Fn(&mut T)>(&mut self, f: F);
    fn fold<F: Fn(T, T) -> T>(&self, init: T, f: F) -> T;
    fn filter<F: Fn(&T) -> bool>(&self, f: F) -> Self;
    fn take(&self, n: usize) -> Self;
    fn skip(&self, n: usize) -> Self;
    fn take_while<F: Fn(&T) -> bool>(&self, f: F) -> Self;
    fn skip_while<F: Fn(&T) -> bool>(&self, f: F) -> Self;
    fn zip_with<F: Fn(T, T) -> T>(&self, f: F, other: &Self) -> Self;
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
    };
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

            fn map<F: Fn($type) -> $type>(&self, f: F) -> Self {
                let v: Vec<$type> = self.to_vec();
                Series::new(v.into_iter().map(f).collect::<Vec<$type>>())
            }

            fn mut_map<F: Fn(&mut $type)>(&mut self, f: F) {
                let v = self.as_slice_mut();
                v.iter_mut().for_each(f);
            }

            fn fold<F: Fn($type, $type) -> $type>(&self, init: $type, f: F) -> $type {
                let v: Vec<$type> = self.to_vec();
                v.into_iter().fold(init, f)
            }

            fn filter<F: Fn(&$type) -> bool>(&self, f: F) -> Self {
                let v: Vec<$type> = self.to_vec();
                Series::new(v.into_iter().filter(|x| f(x)).collect::<Vec<$type>>())
            }

            fn take(&self, n: usize) -> Self {
                let v: Vec<$type> = self.to_vec();
                Series::new(v.into_iter().take(n).collect::<Vec<$type>>())
            }

            fn skip(&self, n: usize) -> Self {
                let v: Vec<$type> = self.to_vec();
                Series::new(v.into_iter().skip(n).collect::<Vec<$type>>())
            }

            fn take_while<F: Fn(&$type) -> bool>(&self, f: F) -> Self {
                let v: Vec<$type> = self.to_vec();
                Series::new(v.into_iter().take_while(|x| f(x)).collect::<Vec<$type>>())
            }

            fn skip_while<F: Fn(&$type) -> bool>(&self, f: F) -> Self {
                let v: Vec<$type> = self.to_vec();
                Series::new(v.into_iter().skip_while(|x| f(x)).collect::<Vec<$type>>())
            }

            fn zip_with<F: Fn($type, $type) -> $type>(&self, f: F, other: &Self) -> Self {
                let v: Vec<$type> = self.to_vec();
                let w: Vec<$type> = other.to_vec();
                Series::new(
                    v.into_iter()
                        .zip(w.into_iter())
                        .map(|(x, y)| f(x, y))
                        .collect::<Vec<$type>>(),
                )
            }
        }
    };
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

    (N; $dtype:expr, $value:expr, $wrapper:expr) => {{
        match $dtype {
            U8 => dtype_case!(u8, $value, $wrapper),
            U16 => dtype_case!(u16, $value, $wrapper),
            U32 => dtype_case!(u32, $value, $wrapper),
            U64 => dtype_case!(u64, $value, $wrapper),
            I8 => dtype_case!(i8, $value, $wrapper),
            I16 => dtype_case!(i16, $value, $wrapper),
            I32 => dtype_case!(i32, $value, $wrapper),
            I64 => dtype_case!(i64, $value, $wrapper),
            F32 => dtype_case!(f32, $value, $wrapper),
            F64 => dtype_case!(f64, $value, $wrapper),
            _ => panic!("Can't use {} to numeric", $dtype);
        }
    }};

    (N; $dtype:expr, $value:expr, $wrapper:expr; $functor:ident) => {{
        match $dtype {
            U8 => dtype_case!($functor<u8>, $value, $wrapper),
            U16 => dtype_case!($functor<u16>, $value, $wrapper),
            U32 => dtype_case!($functor<u32>, $value, $wrapper),
            U64 => dtype_case!($functor<u64>, $value, $wrapper),
            I8 => dtype_case!($functor<i8>, $value, $wrapper),
            I16 => dtype_case!($functor<i16>, $value, $wrapper),
            I32 => dtype_case!($functor<i32>, $value, $wrapper),
            I64 => dtype_case!($functor<i64>, $value, $wrapper),
            F32 => dtype_case!($functor<f32>, $value, $wrapper),
            F64 => dtype_case!($functor<f64>, $value, $wrapper),
            _ => panic!("Can't use {} to numeric", $dtype),
        }
    }};
}

macro_rules! set_space {
    ($elem:expr) => {{
        match $elem.dtype {
            F32 => {
                let elem: f32 = $elem.unwrap();
                let st1 = elem.fmt_lower_exp(2);
                let st2 = elem.to_string();

                if st1.len() < st2.len() {
                    st1
                } else {
                    st2
                }
            }
            F64 => {
                let elem: f64 = $elem.unwrap();
                let st1 = elem.fmt_lower_exp(2);
                let st2 = elem.to_string();

                if st1.len() < st2.len() {
                    st1
                } else {
                    st2
                }
            }
            _ => $elem.to_string(),
        }
    }};

    ($elem:expr, $space:expr) => {{
        match $elem.dtype {
            F32 => {
                let elem: f32 = $elem.unwrap();
                $space = max(
                    $space,
                    min(elem.fmt_lower_exp(2).len(), elem.to_string().len()),
                );
            }
            F64 => {
                let elem: f64 = $elem.unwrap();
                $space = max(
                    $space,
                    min(elem.fmt_lower_exp(2).len(), elem.to_string().len()),
                );
            }
            _ => {
                $space = max($space, $elem.to_string().len());
            }
        }
    }};
}

macro_rules! format_float_vec {
    ($self:expr) => {{
        let mut result = String::new();
        result.push_str("[");
        for i in 0..$self.len() {
            let st1 = $self[i].fmt_lower_exp(2);
            let st2 = $self[i].to_string();
            let st = if st1.len() < st2.len() { st1 } else { st2 };
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
    }};
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
            _ => panic!("Can't convert to {}", $dt2),
        }
    }};
}

macro_rules! dtype_cast_vec {
    ($dt1:expr, $dt2:expr, $to_vec:expr, $wrapper:expr) => {{
        match $dt1 {
            USIZE => dtype_cast_vec_part!(usize, $dt2, $to_vec, $wrapper),
            U8 => match $dt2 {
                Bool => {
                    let y: Vec<u8> = $to_vec;
                    let x: Vec<bool> = y.into_iter().map(|x| x != 0).collect();
                    $wrapper(x)
                }
                Char => {
                    let y: Vec<u8> = $to_vec;
                    let x: Vec<char> = y.into_iter().map(|x| x as char).collect();
                    $wrapper(x)
                }
                _ => dtype_cast_vec_part!(u8, $dt2, $to_vec, $wrapper),
            },
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
            Char => match $dt2 {
                Str => string_cast_vec!(char, $to_vec, $wrapper),
                U8 => {
                    let y: Vec<char> = $to_vec;
                    let x: Vec<u8> = y.into_iter().map(|x| x as u8).collect();
                    $wrapper(x)
                }
                _ => panic!("Can't convert char type to {}", $dt2),
            },
            Bool => match $dt2 {
                U8 => {
                    let y: Vec<bool> = $to_vec;
                    let x: Vec<u8> = y.into_iter().map(|x| if x { 1 } else { 0 }).collect();
                    $wrapper(x)
                }
                Bool => {
                    let y: Vec<bool> = $to_vec;
                    $wrapper(y)
                }
                _ => panic!("Can't convert bool type to {}", $dt2),
            },
        }
    }};
}

fn len<T>(x: Vec<T>) -> usize {
    x.len()
}

fn to_string<T: fmt::Display>(x: T) -> String {
    x.to_string()
}

#[cfg(feature = "nc")]
fn dtype_to_vtype(dt: DType) -> netcdf::types::BasicType {
    match dt {
        USIZE => netcdf::types::BasicType::Uint64,
        U8 => netcdf::types::BasicType::Ubyte,
        U16 => netcdf::types::BasicType::Ushort,
        U32 => netcdf::types::BasicType::Uint,
        U64 => netcdf::types::BasicType::Uint64,
        ISIZE => netcdf::types::BasicType::Int64,
        I8 => netcdf::types::BasicType::Byte,
        I16 => netcdf::types::BasicType::Short,
        I32 => netcdf::types::BasicType::Int,
        I64 => netcdf::types::BasicType::Int64,
        F32 => netcdf::types::BasicType::Float,
        F64 => netcdf::types::BasicType::Double,
        Bool => netcdf::types::BasicType::Ubyte,
        Char => netcdf::types::BasicType::Ubyte,
        _ => panic!("Can't convert type to netcdf::types::BasicType"),
    }
}

#[cfg(feature = "nc")]
fn vtype_to_dtype(dv: netcdf::types::BasicType) -> DType {
    match dv {
        netcdf::types::BasicType::Ubyte => U8,
        netcdf::types::BasicType::Ushort => U16,
        netcdf::types::BasicType::Uint => U32,
        netcdf::types::BasicType::Uint64 => U64,
        netcdf::types::BasicType::Byte => I8,
        netcdf::types::BasicType::Short => I16,
        netcdf::types::BasicType::Int => I32,
        netcdf::types::BasicType::Int64 => I64,
        netcdf::types::BasicType::Float => F32,
        netcdf::types::BasicType::Double => F64,
        netcdf::types::BasicType::Char => Char,
    }
}

#[cfg(feature = "nc")]
fn nc_put_value<T: Numeric>(var: &mut VariableMut, v: Vec<T>) -> Result<(), netcdf::error::Error> {
    var.put_values(&v, None, None)
}

#[cfg(feature = "nc")]
fn nc_read_value<T: Numeric + Default + Clone>(
    val: &Variable,
    v: Vec<T>,
) -> Result<Series, netcdf::error::Error>
where
    Series: TypedVector<T>,
{
    let mut v = v;
    v.resize_with(val.len(), Default::default);
    val.values_to(&mut v, None, None)?;
    Ok(Series::new(v.clone()))
}

#[cfg(feature = "parquet")]
fn dtype_to_arrow(dt: DType) -> DataType {
    match dt {
        USIZE => DataType::UInt64,
        U8 => DataType::UInt8,
        U16 => DataType::UInt16,
        U32 => DataType::UInt32,
        U64 => DataType::UInt64,
        ISIZE => DataType::Int64,
        I8 => DataType::Int8,
        I16 => DataType::Int16,
        I32 => DataType::Int32,
        I64 => DataType::Int64,
        F32 => DataType::Float32,
        F64 => DataType::Float64,
        Bool => DataType::Boolean,
        Str => DataType::Utf8,
        Char => DataType::Utf8,
    }
}

#[cfg(feature = "parquet")]
fn arrow_to_dtype(dt: DataType) -> DType {
    match dt {
        DataType::Boolean => Bool,
        DataType::Int8 => I8,
        DataType::Int16 => I16,
        DataType::Int32 => I32,
        DataType::Int64 => I64,
        DataType::UInt8 => U8,
        DataType::UInt16 => U16,
        DataType::UInt32 => U32,
        DataType::UInt64 => U64,
        // DataType::Float16 => DType::F16,
        DataType::Float32 => F32,
        DataType::Float64 => F64,
        DataType::Utf8 => Str,
        _ => unimplemented!(),
    }
}

#[cfg(feature = "parquet")]
macro_rules! dtype_case_to_arrow {
    ($ty:ty, $to_arr:expr, $value:expr, $chunk_vec:expr; $length:expr) => {{
        let v: Vec<$ty> = $value;
        let v_wrap = (0usize..$length)
            .map(|i| {
                if i < v.len() {
                    Some(v[i].clone())
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        let arr = $to_arr(v_wrap);
        $chunk_vec.push(Arc::from(arr) as Arc<dyn Array>);
    }};
}

#[cfg(feature = "parquet")]
macro_rules! dtype_match_to_arrow {
    ($dtype:expr, $value:expr, $chunk_vec:expr; $length:expr) => {{
        match $dtype {
            Bool => dtype_case_to_arrow!(bool, BooleanArray::from, $value, $chunk_vec; $length),
            Str => dtype_case_to_arrow!(String, StringArray::from, $value, $chunk_vec; $length),
            Char => {
                let v: Vec<char> = $value;
                let v = v.into_iter().map(|t| t.to_string()).collect::<Vec<_>>();
                dtype_case_to_arrow!(String, StringArray::from, v, $chunk_vec; $length)
            }
            USIZE => dtype_case_to_arrow!(u64, PrimitiveArray::<UInt64Type>::from, $value, $chunk_vec; $length),
            U8 => dtype_case_to_arrow!(u8, PrimitiveArray::<UInt8Type>::from, $value, $chunk_vec; $length),
            U16 => dtype_case_to_arrow!(u16, PrimitiveArray::<UInt16Type>::from, $value, $chunk_vec; $length),
            U32 => dtype_case_to_arrow!(u32, PrimitiveArray::<UInt32Type>::from, $value, $chunk_vec; $length),
            U64 => dtype_case_to_arrow!(u64, PrimitiveArray::<UInt64Type>::from, $value, $chunk_vec; $length),
            ISIZE => dtype_case_to_arrow!(i64, PrimitiveArray::<Int64Type>::from, $value, $chunk_vec; $length),
            I8 => dtype_case_to_arrow!(i8, PrimitiveArray::<Int8Type>::from, $value, $chunk_vec; $length),
            I16 => dtype_case_to_arrow!(i16, PrimitiveArray::<Int16Type>::from, $value, $chunk_vec; $length),
            I32 => dtype_case_to_arrow!(i32, PrimitiveArray::<Int32Type>::from, $value, $chunk_vec; $length),
            I64 => dtype_case_to_arrow!(i64, PrimitiveArray::<Int64Type>::from, $value, $chunk_vec; $length),
            F32 => dtype_case_to_arrow!(f32, PrimitiveArray::<Float32Type>::from, $value, $chunk_vec; $length),
            F64 => dtype_case_to_arrow!(f64, PrimitiveArray::<Float64Type>::from, $value, $chunk_vec; $length),
        }
    }};
}

fn add_vec<T: std::ops::Add<T, Output = T> + Clone>(v: Vec<T>, w: Vec<T>) -> Series
where
    Series: TypedVector<T>,
{
    Series::new(v.into_iter().zip(w).map(|(x, y)| x + y).collect::<Vec<T>>())
}

fn sub_vec<T: std::ops::Sub<T, Output = T> + Clone>(v: Vec<T>, w: Vec<T>) -> Series
where
    Series: TypedVector<T>,
{
    Series::new(v.into_iter().zip(w).map(|(x, y)| x - y).collect::<Vec<T>>())
}

fn mul_scalar<T: std::ops::Mul<T, Output = T> + Clone + Copy>(v: Vec<T>, s: T) -> Series
where
    Series: TypedVector<T>,
{
    Series::new(v.into_iter().map(|x| x * s).collect::<Vec<T>>())
}

// =============================================================================
// Implementations of DType variables
// =============================================================================
impl DType {
    /// Check for static numeric type
    pub fn is_numeric(&self) -> bool {
        match self {
            Bool => false,
            Str => false,
            Char => false,
            USIZE => false,
            ISIZE => false,
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
    /// Scalar to length 1 Series
    pub fn to_series(self) -> Series {
        dtype_match!(self.dtype, vec![self.unwrap()], Series::new; Vec)
    }

    pub fn to_string(self) -> String {
        dtype_match!(self.dtype, self.unwrap(), to_string)
    }
}

impl Series {
    /// Getter for Series
    ///
    /// # Examples
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = Series::new(vec![1i32,2,3,4]);
    ///     let x = a.at(0);
    ///
    ///     assert_eq!(x, Scalar::new(1i32));
    /// }
    /// ```
    pub fn at(&self, i: usize) -> Scalar {
        dtype_match!(self.dtype, self.at_raw(i), Scalar::new)
    }

    /// Length for Series
    pub fn len(&self) -> usize {
        dtype_match!(self.dtype, self.as_slice().to_vec(), len; Vec)
    }

    /// Explicit type casting for Series
    pub fn to_type(&self, dtype: DType) -> Series {
        dtype_cast_vec!(self.dtype, dtype, self.to_vec(), Series::new)
    }

    /// Type casting for Series
    ///
    /// # Examples
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let mut a = Series::new(vec![1i32, 2, 3, 4]);
    ///     a.as_type(USIZE);
    ///     
    ///     assert_eq!(a, Series::new(vec![1usize, 2, 3, 4]));
    /// }
    /// ```
    pub fn as_type(&mut self, dtype: DType) {
        let x = self.to_type(dtype);
        self.dtype = x.dtype;
        self.values = x.values;
    }

    /// Select elements by indices, returning a new Series
    ///
    /// # Examples
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = Series::new(vec![10, 20, 30, 40, 50]);
    ///     let b = a.select_indices(&[0, 2, 4]);
    ///     assert_eq!(b, Series::new(vec![10, 30, 50]));
    /// }
    /// ```
    pub fn select_indices(&self, indices: &[usize]) -> Series {
        macro_rules! extract_by_indices {
            ($array:expr, $type:ty) => {{
                let values: Vec<$type> = indices.iter().map(|&i| $array[i].clone()).collect();
                Series::new(values)
            }};
        }

        match &self.values {
            DTypeArray::USIZE(v) => extract_by_indices!(v, usize),
            DTypeArray::U8(v) => extract_by_indices!(v, u8),
            DTypeArray::U16(v) => extract_by_indices!(v, u16),
            DTypeArray::U32(v) => extract_by_indices!(v, u32),
            DTypeArray::U64(v) => extract_by_indices!(v, u64),
            DTypeArray::ISIZE(v) => extract_by_indices!(v, isize),
            DTypeArray::I8(v) => extract_by_indices!(v, i8),
            DTypeArray::I16(v) => extract_by_indices!(v, i16),
            DTypeArray::I32(v) => extract_by_indices!(v, i32),
            DTypeArray::I64(v) => extract_by_indices!(v, i64),
            DTypeArray::F32(v) => extract_by_indices!(v, f32),
            DTypeArray::F64(v) => extract_by_indices!(v, f64),
            DTypeArray::Bool(v) => extract_by_indices!(v, bool),
            DTypeArray::Str(v) => extract_by_indices!(v, String),
            DTypeArray::Char(v) => extract_by_indices!(v, char),
        }
    }

    /// Convert numeric Series to `Vec<f64>`
    ///
    /// Supports all integer and float types. Non-numeric types (Bool, Char, Str) return an error.
    pub fn to_f64_vec(&self) -> anyhow::Result<Vec<f64>> {
        match self.dtype {
            Bool | Char | Str => anyhow::bail!("Cannot convert {} Series to f64", self.dtype),
            _ => {
                let converted = self.to_type(F64);
                Ok(TypedVector::<f64>::to_vec(&converted))
            }
        }
    }

    // =========================================================================
    // Statistics
    // =========================================================================

    /// Sum of all elements (numeric types only)
    pub fn sum(&self) -> anyhow::Result<f64> {
        let v = self.to_f64_vec()?;
        Ok(v.iter().sum())
    }

    /// Mean of all elements (numeric types only, Welford's algorithm)
    pub fn mean(&self) -> anyhow::Result<f64> {
        use crate::statistics::stat::Statistics;
        let v = self.to_f64_vec()?;
        anyhow::ensure!(!v.is_empty(), "Cannot compute mean of empty Series");
        Ok(v.mean())
    }

    /// Variance of all elements (numeric types only, sample variance)
    pub fn var(&self) -> anyhow::Result<f64> {
        use crate::statistics::stat::Statistics;
        let v = self.to_f64_vec()?;
        anyhow::ensure!(v.len() > 1, "Cannot compute variance of Series with fewer than 2 elements");
        Ok(v.var())
    }

    /// Standard deviation of all elements (numeric types only)
    pub fn sd(&self) -> anyhow::Result<f64> {
        use crate::statistics::stat::Statistics;
        let v = self.to_f64_vec()?;
        anyhow::ensure!(v.len() > 1, "Cannot compute sd of Series with fewer than 2 elements");
        Ok(v.sd())
    }

    /// Minimum value, preserving original type
    pub fn min(&self) -> anyhow::Result<Scalar> {
        anyhow::ensure!(self.len() > 0, "Cannot compute min of empty Series");

        macro_rules! typed_min {
            ($v:expr, $dtype:ident) => {{
                let min_val = $v.iter().cloned().reduce(|a, b| if a <= b { a } else { b }).unwrap();
                Ok(Scalar { value: DTypeValue::$dtype(min_val), dtype: DType::$dtype })
            }};
        }

        match &self.values {
            DTypeArray::USIZE(v) => typed_min!(v, USIZE),
            DTypeArray::U8(v) => typed_min!(v, U8),
            DTypeArray::U16(v) => typed_min!(v, U16),
            DTypeArray::U32(v) => typed_min!(v, U32),
            DTypeArray::U64(v) => typed_min!(v, U64),
            DTypeArray::ISIZE(v) => typed_min!(v, ISIZE),
            DTypeArray::I8(v) => typed_min!(v, I8),
            DTypeArray::I16(v) => typed_min!(v, I16),
            DTypeArray::I32(v) => typed_min!(v, I32),
            DTypeArray::I64(v) => typed_min!(v, I64),
            DTypeArray::F32(v) => typed_min!(v, F32),
            DTypeArray::F64(v) => typed_min!(v, F64),
            DTypeArray::Bool(v) => typed_min!(v, Bool),
            DTypeArray::Char(v) => typed_min!(v, Char),
            DTypeArray::Str(v) => typed_min!(v, Str),
        }
    }

    /// Maximum value, preserving original type
    pub fn max(&self) -> anyhow::Result<Scalar> {
        anyhow::ensure!(self.len() > 0, "Cannot compute max of empty Series");

        macro_rules! typed_max {
            ($v:expr, $dtype:ident) => {{
                let max_val = $v.iter().cloned().reduce(|a, b| if a >= b { a } else { b }).unwrap();
                Ok(Scalar { value: DTypeValue::$dtype(max_val), dtype: DType::$dtype })
            }};
        }

        match &self.values {
            DTypeArray::USIZE(v) => typed_max!(v, USIZE),
            DTypeArray::U8(v) => typed_max!(v, U8),
            DTypeArray::U16(v) => typed_max!(v, U16),
            DTypeArray::U32(v) => typed_max!(v, U32),
            DTypeArray::U64(v) => typed_max!(v, U64),
            DTypeArray::ISIZE(v) => typed_max!(v, ISIZE),
            DTypeArray::I8(v) => typed_max!(v, I8),
            DTypeArray::I16(v) => typed_max!(v, I16),
            DTypeArray::I32(v) => typed_max!(v, I32),
            DTypeArray::I64(v) => typed_max!(v, I64),
            DTypeArray::F32(v) => typed_max!(v, F32),
            DTypeArray::F64(v) => typed_max!(v, F64),
            DTypeArray::Bool(v) => typed_max!(v, Bool),
            DTypeArray::Char(v) => typed_max!(v, Char),
            DTypeArray::Str(v) => typed_max!(v, Str),
        }
    }
}

impl Vector for Series {
    type Scalar = Scalar;

    /// Add series
    ///
    /// # Example
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = Series::new(vec![1,2,3]);
    ///     let b = Series::new(vec![3,2,1]);
    ///     let c = a.add_vec(&b);
    ///     assert_eq!(c, Series::new(vec![4,4,4]));
    /// }
    /// ```
    fn add_vec(&self, rhs: &Self) -> Self {
        assert_eq!(self.dtype, rhs.dtype, "DTypes are not same (add_vec)");
        dtype_match!(
            N;
            self.dtype,
            self.to_vec(),
            |x| add_vec(x, rhs.to_vec());
            Vec
        )
    }

    /// Sub series
    ///
    /// # Example
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = Series::new(vec![4,5,6]);
    ///     let b = Series::new(vec![1,2,3]);
    ///     let c = a.sub_vec(&b);
    ///     assert_eq!(c, Series::new(vec![3,3,3]));
    /// }
    /// ```
    fn sub_vec(&self, rhs: &Self) -> Self {
        assert_eq!(self.dtype, rhs.dtype, "DTypes are not same (add_vec)");
        dtype_match!(
            N;
            self.dtype,
            self.to_vec(),
            |x| sub_vec(x, rhs.to_vec());
            Vec
        )
    }

    /// Mul Scalar
    ///
    /// # Example
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = Series::new(vec![1,2,3]);
    ///     let b = Scalar::new(2);
    ///     let c = a.mul_scalar(b);
    ///     assert_eq!(c, Series::new(vec![2,4,6]));
    /// }
    /// ```
    fn mul_scalar(&self, rhs: Self::Scalar) -> Self {
        assert_eq!(self.dtype, rhs.dtype, "DTypes are not same (mul_scalar)");

        dtype_match!(
            N;
            self.dtype,
            self.to_vec(),
            |x| mul_scalar(x, rhs.unwrap());
            Vec
        )
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

// impl FPVector for Series {
//     type Scalar = Scalar;
//
//     fn fmap<F>(&self, f: F) -> Self where
//         F: Fn(Self::Scalar) -> Self::Scalar {
//         dtype_match!(
//             self.dtype,
//             self.to_vec(),
//             |x| map(x, f);
//             Vec
//         )
//     }
//
//     fn reduce<F, T>(&self, init: T, f: F) -> Self::Scalar where
//         F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar,
//         T: Into<Self::Scalar> {
//         dtype_match!(
//             self.dtype,
//             self.to_vec(),
//             |x| reduce(x, f);
//             Vec
//         )
//     }
//
//     fn zip_with<F>(&self, f: F, other: &Self) -> Self where
//         F: Fn(Self::Scalar, Self::Scalar) -> Self::Scalar {
//         dtype_match!(
//             self.dtype,
//             self.to_vec(),
//             |x| zip_with(x, other.to_vec(), f);
//             Vec
//         )
//     }
//
//     fn filter<F>(&self, f: F) -> Self where
//         F: Fn(Self::Scalar) -> bool {
//         dtype_match!(
//             self.dtype,
//             self.to_vec(),
//             |x| filter(x, f);
//             Vec
//         )
//     }
//
//     fn take(&self, n: usize) -> Self {
//         dtype_match!(
//             self.dtype,
//             self.to_vec(),
//             |x| take(x, n);
//             Vec
//         )
//     }
//
//     fn skip(&self, n: usize) -> Self {
//         dtype_match!(
//             self.dtype,
//             self.to_vec(),
//             |x| skip(x, n);
//             Vec
//         )
//     }
//
//     fn sum(&self) -> Self::Scalar {
//         todo!()
//     }
//
//     fn prod(&self) -> Self::Scalar {
//         todo!()
//     }
// }

// =============================================================================
// Implementation for DataFrame
// =============================================================================

impl DataFrame {
    /// Declare new DataFrame with `Vec<Series>`
    pub fn new(v: Vec<Series>) -> Self {
        let ics = (0usize..v.len()).map(|x| x.to_string()).collect();

        Self { data: v, ics }
    }

    pub fn header(&self) -> &Vec<String> {
        &self.ics
    }

    pub fn header_mut(&mut self) -> &mut Vec<String> {
        &mut self.ics
    }

    /// Change header
    pub fn set_header(&mut self, new_header: Vec<&str>) {
        assert_eq!(self.ics.len(), new_header.len(), "Improper Header length!");
        self.ics = new_header.into_iter().map(|x| x.to_string()).collect();
    }

    /// Push new pair of head, Series to DataFrame
    pub fn push(&mut self, name: &str, series: Series) {
        if !self.ics.is_empty() {
            assert_eq!(
                self.ics.iter().find(|x| x.as_str() == name),
                None,
                "Repetitive index!"
            );
        }
        self.ics.push(name.to_string());
        self.data.push(series);
    }

    /// Extract specific row as DataFrame
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
        let r: usize = self
            .data
            .iter()
            .fold(0, |max_len, column| max(max_len, column.len()));
        let h = self.header();

        let mut result = String::new();

        if r > 100 {
            let lc1 = ((r as f64).log10() as usize) + 5;
            result.push_str(&tab("", lc1));

            let mut space_vec: Vec<usize> = vec![];
            for i in 0..self.data.len() {
                let v = &self[i];
                let mut space = 0usize;
                for j in 0..v.len().min(5) {
                    let elem = v.at(j);
                    set_space!(elem, space);
                }
                if v.len() >= r - 5 {
                    for j in v.len() - 5..v.len() {
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

            for i in 0..5 {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for j in 0..self.data.len() {
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
                result.push('\n');
            }
            result.push_str(&tab("...", lc1));
            for &space in space_vec.iter() {
                result.push_str(&tab("...", space));
            }
            result.push('\n');
            for i in r - 5..r {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for j in 0..self.data.len() {
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
                if i == r - 1 {
                    break;
                }
                result.push('\n');
            }
            return result;
        }

        result.push_str(&tab("", 5));
        let mut space_vec: Vec<usize> = vec![];

        for i in 0..self.data.len() {
            let v = &self[i];
            let mut space = 0usize;
            for j in 0..v.len() {
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

        for i in 0..r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for j in 0..self.data.len() {
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

    /// Type casting for DataFrame
    ///
    /// # Examples
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = Series::new(vec![1i32, 2, 3, 4]);
    ///     let b = Series::new(vec![true, false, false, true]);
    ///     
    ///     let mut df = DataFrame::new(vec![a, b]);    // I32, Bool
    ///     df.as_types(vec![USIZE, U8]);               // USIZE, U8
    ///
    ///     let c = Series::new(vec![1usize, 2, 3, 4]);
    ///     let d = Series::new(vec![1u8, 0, 0, 1]);
    ///     let dg = DataFrame::new(vec![c, d]);
    ///
    ///     assert_eq!(df, dg);
    /// }
    /// ```
    pub fn as_types(&mut self, dtypes: Vec<DType>) {
        assert_eq!(
            self.data.len(),
            dtypes.len(),
            "Length of dtypes are not compatible with DataFrame"
        );
        for (i, dtype) in dtypes.into_iter().enumerate() {
            self[i].as_type(dtype);
        }
    }

    /// Drop specific column by header
    ///
    /// # Examples
    ///
    /// ```rust
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let a = Series::new(vec![1,2,3,4]);
    ///     let b = Series::new(vec![5,6,7,8]);
    ///
    ///     let mut df = DataFrame::new(vec![a.clone(), b]);
    ///     df.set_header(vec!["a", "b"]);
    ///
    ///     let mut dg = DataFrame::new(vec![a]);
    ///     dg.set_header(vec!["a"]);
    ///
    ///     df.drop("b");
    ///
    ///     assert_eq!(df, dg);
    /// }
    /// ```
    pub fn drop(&mut self, col_header: &str) {
        match self.ics.iter().position(|h| h == col_header) {
            Some(index) => {
                self.data.remove(index);
                self.ics.remove(index);
            }
            None => panic!("Can't drop header '{}'", col_header),
        }
    }

    /// Filter DataFrame by specific column
    pub fn filter_by<F>(&self, column: &str, predicate: F) -> anyhow::Result<DataFrame>
    where
        F: Fn(Scalar) -> bool,
    {
        let series = match self.ics.iter().position(|x| x.as_str() == column) {
            Some(i) => &self.data[i],
            None => anyhow::bail!("Column '{}' not found in DataFrame", column),
        };

        let mut indices = Vec::new();
        for i in 0..series.len() {
            let value = series.at(i);
            if predicate(value) {
                indices.push(i);
            }
        }

        let mut new_df = DataFrame::new(vec![]);
        for (col_idx, col_series) in self.data.iter().enumerate() {
            let filtered_series = col_series.select_indices(&indices);
            new_df.push(&self.ics[col_idx], filtered_series);
        }

        Ok(new_df)
    }

    /// Mask DataFrame with a boolean Series
    pub fn mask(&self, mask: &Series) -> anyhow::Result<DataFrame> {
        if mask.len() != self.data[0].len() {
            anyhow::bail!(
                "Mask length ({}) does not match DataFrame row count ({})",
                mask.len(),
                self.data[0].len()
            );
        }

        if mask.dtype != DType::Bool {
            anyhow::bail!("Mask Series must be of type Bool, but got {}", mask.dtype);
        }

        let bool_mask: &[bool] = mask.as_slice();
        let ics: Vec<usize> = bool_mask
            .iter()
            .enumerate()
            .filter_map(|(i, &b)| if b { Some(i) } else { None })
            .collect();

        Ok(self.select_rows(&ics))
    }

    /// Select rows based on indices
    pub fn select_rows(&self, indices: &[usize]) -> DataFrame {
        let mut new_df = DataFrame::new(vec![]);
        for (col_idx, col_series) in self.data.iter().enumerate() {
            let filtered_series = col_series.select_indices(indices);
            new_df.push(&self.ics[col_idx], filtered_series);
        }
        new_df
    }

    // =========================================================================
    // Shape / Info
    // =========================================================================

    /// Number of rows (max column length)
    pub fn nrow(&self) -> usize {
        self.data.iter().fold(0, |acc, s| max(acc, s.len()))
    }

    /// Number of columns
    pub fn ncol(&self) -> usize {
        self.data.len()
    }

    /// Shape as (nrow, ncol)
    pub fn shape(&self) -> (usize, usize) {
        (self.nrow(), self.ncol())
    }

    /// DType of each column
    pub fn dtypes(&self) -> Vec<DType> {
        self.data.iter().map(|s| s.dtype).collect()
    }

    /// Check if the DataFrame has no columns or no rows
    pub fn is_empty(&self) -> bool {
        self.data.is_empty() || self.nrow() == 0
    }

    /// Check if the DataFrame contains a column with the given header
    pub fn contains(&self, col_header: &str) -> bool {
        self.ics.iter().any(|x| x.as_str() == col_header)
    }

    // =========================================================================
    // Row Operations
    // =========================================================================

    /// Return the first `n` rows
    pub fn head(&self, n: usize) -> DataFrame {
        let nrow = self.nrow();
        let end = n.min(nrow);
        let indices: Vec<usize> = (0..end).collect();
        self.select_rows(&indices)
    }

    /// Return the last `n` rows
    pub fn tail(&self, n: usize) -> DataFrame {
        let nrow = self.nrow();
        let start = nrow.saturating_sub(n);
        let indices: Vec<usize> = (start..nrow).collect();
        self.select_rows(&indices)
    }

    /// Return a slice of rows starting at `offset` with the given `length`
    pub fn slice(&self, offset: usize, length: usize) -> DataFrame {
        let nrow = self.nrow();
        let end = (offset + length).min(nrow);
        let indices: Vec<usize> = (offset..end).collect();
        self.select_rows(&indices)
    }

    // =========================================================================
    // Column Operations
    // =========================================================================

    /// Select specific columns by name, returning a new DataFrame
    ///
    /// Panics if any column name does not exist.
    pub fn select(&self, columns: &[&str]) -> DataFrame {
        let mut new_df = DataFrame::new(vec![]);
        for &col in columns {
            let i = self
                .ics
                .iter()
                .position(|x| x.as_str() == col)
                .unwrap_or_else(|| panic!("Column '{}' not found in DataFrame", col));
            new_df.push(col, self.data[i].clone());
        }
        new_df
    }

    /// Rename a column in-place
    ///
    /// Panics if the old column name does not exist.
    pub fn rename(&mut self, old: &str, new: &str) {
        let i = self
            .ics
            .iter()
            .position(|x| x.as_str() == old)
            .unwrap_or_else(|| panic!("Column '{}' not found in DataFrame", old));
        self.ics[i] = new.to_string();
    }

    /// Return column names as `Vec<&str>`
    pub fn column_names(&self) -> Vec<&str> {
        self.ics.iter().map(|s| s.as_str()).collect()
    }

    /// Select columns whose dtype is in the given list
    pub fn select_dtypes(&self, dtypes: &[DType]) -> DataFrame {
        let mut new_df = DataFrame::new(vec![]);
        for (i, series) in self.data.iter().enumerate() {
            if dtypes.contains(&series.dtype) {
                new_df.push(&self.ics[i], series.clone());
            }
        }
        new_df
    }

    // =========================================================================
    // DataFrame-level Statistics
    // =========================================================================

    /// Compute descriptive statistics for numeric columns
    ///
    /// Returns a DataFrame with rows: count, mean, sd, min, max
    /// and one column per numeric column from the original DataFrame.
    pub fn describe(&self) -> DataFrame {
        use crate::statistics::stat::Statistics;

        let stat_labels = vec!["count", "mean", "sd", "min", "max"];
        let mut result = DataFrame::new(vec![]);
        result.push("stat", Series::new(stat_labels.iter().map(|s| s.to_string()).collect::<Vec<String>>()));

        for (i, series) in self.data.iter().enumerate() {
            if let Ok(v) = series.to_f64_vec() {
                if v.is_empty() {
                    continue;
                }
                let count = v.len() as f64;
                let mean = v.mean();
                let sd = if v.len() > 1 { v.sd() } else { 0.0 };
                let min_val = v.iter().cloned().fold(f64::INFINITY, f64::min);
                let max_val = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                result.push(
                    &self.ics[i],
                    Series::new(vec![count, mean, sd, min_val, max_val]),
                );
            }
        }

        result
    }

    /// Sum of each numeric column as a single-row DataFrame
    pub fn sum(&self) -> DataFrame {
        let mut result = DataFrame::new(vec![]);
        for (i, series) in self.data.iter().enumerate() {
            if let Ok(v) = series.to_f64_vec() {
                let s: f64 = v.iter().sum();
                result.push(&self.ics[i], Series::new(vec![s]));
            }
        }
        result
    }

    /// Mean of each numeric column as a single-row DataFrame
    pub fn mean(&self) -> DataFrame {
        use crate::statistics::stat::Statistics;

        let mut result = DataFrame::new(vec![]);
        for (i, series) in self.data.iter().enumerate() {
            if let Ok(v) = series.to_f64_vec() {
                if v.is_empty() {
                    continue;
                }
                let m = v.mean();
                result.push(&self.ics[i], Series::new(vec![m]));
            }
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

// =============================================================================
// IO Implementations
// =============================================================================

/// To handle CSV file format
#[cfg(feature = "csv")]
pub trait WithCSV: Sized {
    fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>>;
    fn read_csv(file_path: &str, delimiter: char) -> Result<Self, Box<dyn Error>>;
}

#[cfg(feature = "csv")]
impl WithCSV for DataFrame {
    /// Write csv file
    fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r: usize = self
            .data
            .iter()
            .fold(0, |max_len, column| max(max_len, column.len()));
        let c: usize = self.data.len();
        wtr.write_record(self.header().clone())?;

        for i in 0..r {
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

    /// Read csv file with delimiter
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

/// To handle with NetCDF file format
#[cfg(feature = "nc")]
pub trait WithNetCDF: Sized {
    fn write_nc(&self, file_path: &str) -> Result<(), Box<dyn Error>>;
    fn read_nc(file_path: &str) -> Result<Self, Box<dyn Error>>;
    fn read_nc_by_header(file_path: &str, header: Vec<&str>) -> Result<Self, Box<dyn Error>>;
}

#[cfg(feature = "nc")]
impl WithNetCDF for DataFrame {
    /// write netcdf file
    fn write_nc(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut f = netcdf::create(file_path)?;

        for (i, h) in self.header().iter().enumerate() {
            let dim_name = format!("{}th col", i);
            let v = &self[h.as_str()];
            let dim = v.len();
            f.add_dimension(&dim_name, dim)?;
            match v.dtype {
                dtype if dtype.is_numeric() => {
                    let vtype = dtype_to_vtype(dtype);
                    let var = &mut f.add_variable_with_type(
                        h,
                        &[&dim_name],
                        &VariableType::Basic(vtype),
                    )?;
                    dtype_match!(N; dtype, v.to_vec(), |v| nc_put_value(var, v); Vec)?;
                }
                Str => {
                    let var = &mut f.add_string_variable(h, &[&dim_name])?;
                    let v_s: &[String] = v.as_slice();
                    for (i, s) in v_s.iter().enumerate() {
                        var.put_string(s, Some(&[i]))?;
                    }
                }
                USIZE => {
                    let v = v.to_type(U64);
                    let var = &mut f.add_variable::<u64>(h, &[&dim_name])?;
                    let v_slice: &[u64] = v.as_slice();
                    var.put_values(v_slice, None, None)?;
                }
                ISIZE => {
                    let v = v.to_type(I64);
                    let var = &mut f.add_variable::<i64>(h, &[&dim_name])?;
                    let v_slice: &[i64] = v.as_slice();
                    var.put_values(v_slice, None, None)?;
                }
                Bool => {
                    let v = v.to_type(U8);
                    let var = &mut f.add_variable::<u8>(h, &[&dim_name])?;
                    let v_slice: &[u8] = v.as_slice();
                    var.put_values(v_slice, None, None)?;
                }
                Char => {
                    let v = v.to_type(U8);
                    let var = &mut f.add_variable::<u8>(h, &[&dim_name])?;
                    let v_slice: &[u8] = v.as_slice();
                    var.put_values(v_slice, None, None)?;
                }
                _ => unreachable!(),
            }
        }

        Ok(())
    }

    /// Read netcdf to DataFrame
    fn read_nc(file_path: &str) -> Result<Self, Box<dyn Error>> {
        let f = netcdf::open(file_path)?;
        let mut df = DataFrame::new(vec![]);
        for v in f.variables() {
            let h = v.name();
            if v.vartype().is_string() {
                let mut data: Vec<String> = vec![Default::default(); v.len()];
                for i in 0..v.len() {
                    data[i] = v.string_value(Some(&[i]))?;
                }
                df.push(&h, Series::new(data));
            } else {
                let dtype = vtype_to_dtype(v.vartype().as_basic().unwrap());
                let series = dtype_match!(N; dtype, vec![], |vec| nc_read_value(&v, vec); Vec)?;
                df.push(&h, series);
            }
        }
        Ok(df)
    }

    /// Read netcdf to DataFrame with specific header
    ///
    /// # Example
    ///
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() -> Result<(), Box<dyn Error>> {
    ///     let mut df = DataFrame::new(vec![]);
    ///     df.push("a", Series::new(vec![1,2,3,4]));
    ///     df.push("b", Series::new(vec!['a', 'b', 'c', 'd']));
    ///     df.push("c", Series::new(c!(0.1, 0.2, 0.3, 0.4)));
    ///     df.write_nc("example_data/doc_nc2.nc")?;
    ///
    ///     let dg = DataFrame::read_nc_by_header("example_data/doc_nc2.nc", vec!["a", "c"])?;
    ///
    ///     df.drop("b");
    ///
    ///     assert_eq!(df, dg);
    ///     
    ///     Ok(())
    /// }
    /// ```
    fn read_nc_by_header(file_path: &str, header: Vec<&str>) -> Result<Self, Box<dyn Error>> {
        let f = netcdf::open(file_path)?;
        let mut df = DataFrame::new(vec![]);
        for h in header {
            let v = match f.variable(h) {
                Some(val) => val,
                None => panic!("There are no corresponding values"),
            };
            if v.vartype().is_string() {
                let mut data: Vec<String> = vec![Default::default(); v.len()];
                for i in 0..v.len() {
                    data[i] = v.string_value(Some(&[i]))?;
                }
                df.push(&h, Series::new(data));
            } else {
                let dtype = vtype_to_dtype(v.vartype().as_basic().unwrap());
                let series = dtype_match!(N; dtype, vec![], |vec| nc_read_value(&v, vec); Vec)?;
                df.push(&h, series);
            }
        }
        Ok(df)
    }
}

/// To handle parquet format
#[cfg(feature = "parquet")]
pub trait WithParquet {
    fn write_parquet(
        &self,
        file_path: &str,
        compression: Compression,
    ) -> Result<(), Box<dyn Error>>;
    fn read_parquet(file_path: &str) -> Result<Self, Box<dyn Error>>
    where
        Self: Sized;
    // fn read_parquet_by_header(file_path: &str, header: Vec<&str>) -> Result<Self, Box<dyn Error>> where Self: Sized;
}

/// This macro handles the repetitive logic of processing a column from an Arrow array,
/// converting it to a `Vec<T>`, and then inserting or updating it in the provided HashMap.
///
/// # Arguments
/// - `$hash_map`: The mutable HashMap storing the column series.
/// - `$h`: The column name (header).
/// - `$arr`: The `ArrayRef` (the raw column data from Arrow).
/// - `$arrow_type`: The concrete Arrow array type to downcast to (e.g., `BooleanArray`).
/// - `$rust_type`: The target Rust type for the `Vec` (e.g., `bool`).
/// - `|$concrete_array:ident| $extract_body:expr`: A closure-like expression that defines
///   how to extract the data from the downcasted array into a `Vec<$rust_type>`.
#[cfg(feature = "parquet")]
macro_rules! process_column {
    ($hash_map:expr, $h:expr, $arr:expr, $arrow_type:ty, $rust_type:ty, |$concrete_array:ident| $extract_body:expr) => {{
        // Downcast the generic array to the specific Arrow array type.
        let $concrete_array = $arr.as_any().downcast_ref::<$arrow_type>().unwrap();
        // Apply the provided logic to extract data into a Vec.
        let data: Vec<$rust_type> = $extract_body;

        // Check if the column already exists in the map.
        if let Some(existing_data) = $hash_map.get_mut($h) {
            // If it exists, extend the existing vector with the new data.
            let mut vec_data: Vec<$rust_type> = existing_data.to_vec();
            vec_data.extend(data.iter().cloned());
            $hash_map.insert($h.clone(), Series::new(vec_data));
        } else {
            // If it's a new column, insert a new Series.
            $hash_map.insert($h.clone(), Series::new(data));
        }
    }};
}

#[cfg(feature = "parquet")]
impl WithParquet for DataFrame {
    /// Write DataFrame to parquet
    fn write_parquet(
        &self,
        file_path: &str,
        compression: Compression,
    ) -> Result<(), Box<dyn Error>> {
        let mut schema_vec = vec![];
        let mut arr_vec = vec![];

        let max_length = self.data.iter().fold(0usize, |acc, x| acc.max(x.len()));

        for h in self.header().iter() {
            let v = &self[h.as_str()];
            let field = Field::new(h.as_str(), dtype_to_arrow(v.dtype), false);

            dtype_match_to_arrow!(v.dtype, v.to_vec(), arr_vec; max_length);
            schema_vec.push(field);
        }

        let schema = Arc::new(Schema::new(schema_vec));
        let parquet_schema = ArrowSchemaConverter::new()
            .convert(&schema)
            .map_err(|e| format!("Failed to convert schema: {}", e))?;
        let writer_properties = WriterProperties::builder()
            .set_compression(compression)
            .build();
        let props = Arc::new(writer_properties);

        let col_writers = get_column_writers(&parquet_schema, &props, &schema)?;
        let mut workers: Vec<_> = col_writers
            .into_iter()
            .map(|mut col_writer| {
                let (send, recv) = std::sync::mpsc::channel::<ArrowLeafColumn>();
                let handle = std::thread::spawn(move || {
                    for col in recv {
                        col_writer.write(&col)?;
                    }
                    col_writer.close()
                });
                (handle, send)
            })
            .collect();

        let root_schema = parquet_schema.root_schema_ptr();
        let mut output_file = std::fs::File::create(file_path)?;
        let mut writer = SerializedFileWriter::new(&mut output_file, root_schema, props.clone())?;

        let mut row_group_writer: SerializedRowGroupWriter<'_, _> = writer.next_row_group()?;

        let mut worker_iter = workers.iter_mut();
        for (arr, field) in arr_vec.iter().zip(&schema.fields) {
            for leaves in compute_leaves(field, &Arc::new(arr))? {
                worker_iter.next().unwrap().1.send(leaves)?;
            }
        }

        for (handle, send) in workers {
            use parquet::arrow::arrow_writer::ArrowColumnChunk;

            drop(send);
            let chunk: ArrowColumnChunk = handle.join().unwrap().unwrap();
            chunk.append_to_row_group(&mut row_group_writer)?;
        }
        row_group_writer.close()?;
        writer.close()?;

        Ok(())
    }

    /// Read parquet to DataFrame
    fn read_parquet(file_path: &str) -> Result<Self, Box<dyn Error>>
    where
        Self: Sized,
    {
        use parquet::arrow::arrow_reader::ParquetRecordBatchReader;

        let mut df = DataFrame::new(vec![]);

        let file = std::fs::File::open(file_path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file.try_clone()?)?;
        let schema = builder.schema();
        let fields = schema.fields.clone();
        let mut batch_size = usize::MAX; // Use maximum batch size
        let reader: ParquetRecordBatchReader = loop {
            let builder = ParquetRecordBatchReaderBuilder::try_new(file.try_clone()?)?;
            let reader = builder.with_batch_size(batch_size).build();
            match reader {
                Ok(r) => break r,
                Err(e) => {
                    if batch_size > 0 {
                        batch_size /= 10; // Reduce batch size if error occurs
                    } else {
                        println!(
                            "Failed to read parquet file: {} with eventually batch size 1",
                            e
                        );
                        return Err(Box::new(e));
                    }
                }
            }
        };
        let all_batches: Vec<_> = reader.collect::<Result<Vec<_>, _>>()?;

        let mut hash_map = IndexMap::<String, Series>::new();
        for batch in all_batches {
            let arrs = batch.columns();

            for (field, arr) in fields.iter().zip(arrs) {
                let h = field.name();
                let dt = field.data_type();
                let at = arrow_to_dtype(dt.clone());
                match at {
                    Bool => process_column!(hash_map, h, arr, BooleanArray, bool, |d| d
                        .values()
                        .iter()
                        .collect()),
                    Char => process_column!(hash_map, h, arr, StringArray, char, |d| d
                        .iter()
                        .filter_map(|opt_s| opt_s.and_then(|s| s.chars().next()))
                        .collect()),
                    Str => process_column!(hash_map, h, arr, StringArray, String, |d| d
                        .iter()
                        .filter_map(|opt_s| opt_s.map(String::from))
                        .collect()),
                    USIZE => {
                        process_column!(hash_map, h, arr, PrimitiveArray<UInt64Type>, usize, |d| d
                            .values()
                            .iter()
                            .map(|&x| x as usize)
                            .collect())
                    }
                    U8 => process_column!(hash_map, h, arr, PrimitiveArray<UInt8Type>, u8, |d| d
                        .values()
                        .to_vec()),
                    U16 => {
                        process_column!(hash_map, h, arr, PrimitiveArray<UInt16Type>, u16, |d| d
                            .values()
                            .to_vec())
                    }
                    U32 => {
                        process_column!(hash_map, h, arr, PrimitiveArray<UInt32Type>, u32, |d| d
                            .values()
                            .to_vec())
                    }
                    U64 => {
                        process_column!(hash_map, h, arr, PrimitiveArray<UInt64Type>, u64, |d| d
                            .values()
                            .to_vec())
                    }
                    ISIZE => {
                        process_column!(hash_map, h, arr, PrimitiveArray<Int64Type>, isize, |d| d
                            .values()
                            .iter()
                            .map(|&x| x as isize)
                            .collect())
                    }
                    I8 => process_column!(hash_map, h, arr, PrimitiveArray<Int8Type>, i8, |d| d
                        .values()
                        .to_vec()),
                    I16 => process_column!(hash_map, h, arr, PrimitiveArray<Int16Type>, i16, |d| d
                        .values()
                        .to_vec()),
                    I32 => process_column!(hash_map, h, arr, PrimitiveArray<Int32Type>, i32, |d| d
                        .values()
                        .to_vec()),
                    I64 => process_column!(hash_map, h, arr, PrimitiveArray<Int64Type>, i64, |d| d
                        .values()
                        .to_vec()),
                    F32 => {
                        process_column!(hash_map, h, arr, PrimitiveArray<Float32Type>, f32, |d| d
                            .values()
                            .to_vec())
                    }
                    F64 => {
                        process_column!(hash_map, h, arr, PrimitiveArray<Float64Type>, f64, |d| d
                            .values()
                            .to_vec())
                    }
                }
            }
        }

        for (h, data) in hash_map {
            df.push(&h, data);
        }

        Ok(df)
    }
}
