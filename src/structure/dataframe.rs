//! R-like dataframe
//!
//! ## Declare DataFrame
//!
//! * There are various ways to declare `DataFrame`
//!     * `DataFrame::new()` - Empty dataframe
//!     * `DataFrame::with_header(header: Vec<&str>)` - Empty dataframe with header
//!     * `DataFrame::from_matrix(mat: Matrix)` - from matrix
//!
//! ### Example
//!
//!
//! #### 1. Empty DataFrame
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // Empty DataFrame
//!     let mut a = DataFrame::new();
//!
//!     // Insert method
//!     a.insert("x", vec![1f64, 2f64, 3f64]);
//!     a.insert("y", vec![4f64, 5f64]);
//!     a.insert("z", vec![6f64]);
//!
//!     // Print
//!     a.print();
//!     //         x    y    z
//!     // r[0]    1    4    6
//!     // r[1]    2    5
//!     // r[2]    3
//! }
//! ```
//!
//! #### 2. With Header
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     // With header
//!     let mut a = DataFrame::with_header(vec!["x", "y", "z"]);
//!
//!     // IndexMut method
//!     a["x"] = vec![1f64, 2f64, 3f64];
//!     a["y"] = vec![4f64, 5f64];
//!     a["z"] = vec![6f64];
//!
//!     a.print(); // Same result
//! }
//! ```
//!
//! ## Get Data
//!
//! * `get(&self, key: &str) -> &Vec<f64>`
//! * `Index<&str>`
//!
//! ### Example
//!
//! ```rust
//! #[macro_use]
//! extern crate peroxide;
//! use peroxide::fuga::*;
//!
//! fn main() {
//!     let mut a = DataFrame::with_header(vec!["x", "y", "z"]);
//!     a["x"] = c!(1,2,3);
//!     a["y"] = c!(4,5);
//!     a["z"] = c!(6);
//!
//!     // Get
//!     assert_eq!(&c!(1,2,3), a.get("x").unwrap());
//!     assert_eq!(c!(4,5), a["y"]);
//! }
//! ```
//!
//! ## Read & Write with `netcdf`, `csv` & `json`
//!
//! * Effectiveness : `netcdf` >> `json` > `csv`
//! * File size : `netcdf` > `json` >> `csv`
//!
//! | Command | Mean [s] | Min [s] | Max [s] | Relative |
//! |:---|---:|---:|---:|---:|
//! | `./target/release/df_bench_csv --features dataframe` | 1.027 ± 0.016 | 0.998 | 1.056 | 27.0 |
//! | `./target/release/df_bench_nc --features dataframe` | 0.038 ± 0.002 | 0.035 | 0.042 | 1.0 |
//! | `./target/release/df_bench_json --features dataframe` | 0.652 ± 0.021 | 0.619 | 0.619 | 1.0 |
//!
//! * Benchmark codes are as follows.
//!
//! ### Example (Benchmark codes)
//!
//! #### 1. netcdf
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//! use std::error::Error;
//! use std::fs;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     let mut df = DataFrame::with_header(vec!["x", "y", "z"]);
//!     df["x"] = vec![0f64; 1000_000];
//!     df["y"] = vec![0f64; 1000_000];
//!     df["z"] = vec![0f64; 1000_000];
//!
//!     df.write_nc("example_data/df_bench.nc")?;
//!
//!     // read_nc(path: &str, header: Vec<&str>)
//!     let dg = DataFrame::read_nc("example_data/df_bench.nc")?;
//!     // or read_nc with specific header
//!     let dg = DataFrame::read_nc_by_header("example_data/df_bench.nc", vec!["x", "z"])?;
//!     dg.print();
//!
//!     // For test, we remove nc file
//!     fs::remove_file("example_data/df_bench.nc")?;
//!
//!     Ok(())
//! }
//! ```
//!
//! #### 2. CSV
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//! use std::error::Error;
//! use std::fs; // To remove data file
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     let mut df = DataFrame::with_header(vec!["x", "y", "z"]);
//!     df["x"] = vec![0f64; 1000_000];
//!     df["y"] = vec![0f64; 1000_000];
//!     df["z"] = vec![0f64; 1000_000];
//!
//!     df.write_csv("example_data/df_bench.csv")?;
//!
//!     // read_csv(path: &str, delimiter: char)
//!     let dg = DataFrame::read_csv("example_data/df_bench.csv", ',')?;
//!     dg.print();
//!
//!     // For test, remove csv file
//!     fs::remove_file("example_data/df_bench.csv")?;
//!
//!     Ok(())
//! }
//! ```
//!
//! #### 3. JSON
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::fuga::*;
//! use std::fs::File;
//! use std::error::Error;
//! use std::io::BufWriter;
//! use std::fs;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     let mut df = DataFrame::with_header(vec!["x", "y", "z"]);
//!     df["x"] = vec![0f64; 1000_000];
//!     df["y"] = vec![0f64; 1000_000];
//!     df["z"] = vec![0f64; 1000_000];
//!
//!     let json = df.to_json_value();
//!
//!     // To write JSON
//!     let file = File::create("example_data/df_bench.json")?;
//!     let mut writer = BufWriter::new(file);
//!     json.write(&mut writer)?;
//!
//!     // To read JSON (unimplemented)
//!
//!     // For Test, remove json file
//!     fs::remove_file("example_data/df_bench.json")?;
//!
//!     Ok(())
//! }
//! ```

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
pub struct DataFrame {
    pub data: IndexMap<String, Vec<f64>>,
}

impl PartialEq for DataFrame {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

/// Pretty view for DataFrame
impl fmt::Display for DataFrame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

impl Index<&str> for DataFrame {
    type Output = Vec<f64>;

    fn index(&self, index: &str) -> &Self::Output {
        self.get(index).unwrap()
    }
}

impl IndexMut<&str> for DataFrame {
    fn index_mut(&mut self, index: &str) -> &mut Self::Output {
        match self.data.get_mut(index) {
            Some(v) => v,
            None => panic!("There are no corresponding value"),
        }
    }
}

impl Into<Matrix> for DataFrame {
    fn into(self) -> Matrix {
        self.to_matrix()
    }
}

impl Into<Matrix> for &DataFrame {
    fn into(self) -> Matrix {
        self.to_matrix()
    }
}

#[allow(unused_parens)]
impl DataFrame {
    /// Declare empty DataFrame
    pub fn new() -> Self {
        DataFrame {
            data: IndexMap::new(),
        }
    }

    /// Declare empty DataFrame with header
    pub fn with_header(header: Vec<&str>) -> Self {
        let l = header.len();
        let mut data = IndexMap::with_capacity(l);
        for i in 0..l {
            data.insert(header[i].to_string(), vec![]);
        }
        DataFrame { data }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Insert key & value pair (or only value)
    pub fn insert(&mut self, key: &str, value: Vec<f64>) {
        self.data.insert(key.to_owned(), value);
    }

    /// Insert row
    pub fn insert_row(&mut self, value: Vec<f64>) {
        assert_eq!(self.data.len(), value.len());
        for (v, val) in self.data.values_mut().zip(value) {
            v.push(val);
        }
    }

    /// Get value by ref
    pub fn get(&self, head: &str) -> Option<&Vec<f64>> {
        self.data.get(&head.to_string())
    }

    /// Get iterator of headers
    pub fn headers(&self) -> Keys<String, Vec<f64>> {
        self.data.keys()
    }

    /// Change header
    ///
    /// ```
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// let a = ml_matrix("1 2;3 4");
    /// let mut df = DataFrame::from_matrix(a);
    /// df.print();
    ///
    /// df.set_header(vec!["x", "y"]);
    /// df.print();
    /// ```
    pub fn set_header(&mut self, header: Vec<&str>) {
        let mut im = IndexMap::with_capacity(header.len());
        for ((_, v), &head) in self.data.drain(..).zip(header.iter()) {
            im.insert(head.to_owned(), v);
        }
        self.data = im;
        //for i in 0..header.len() {
        //    match self.data.get_index_mut(i) {
        //        Some((k, _)) => {
        //            *k = header[i].to_string();
        //        }
        //        None => panic!("New header is longer than original header"),
        //    }
        //}
    }

    /// Convert to matrix
    pub fn to_matrix(&self) -> Matrix {
        let mut data: Vec<f64> = vec![];
        let mut r = 0usize;
        let mut c = 0usize;
        self.data.values().for_each(|v| {
            if r == 0 {
                r = v.len();
            } else {
                assert_eq!(r, v.len());
            }
            c += 1;
            data.extend(v);
        });
        matrix(data, r, c, Col)
    }

    /// Convert from matrix
    pub fn from_matrix(mat: Matrix) -> Self {
        let mut df: DataFrame = DataFrame::new();
        for i in 0..mat.col {
            df.insert(format!("{}", i).as_str(), mat.col(i));
        }
        df
    }

    /// For pretty print
    pub fn spread(&self) -> String {
        let r: usize = self
            .data
            .values()
            .fold(0, |max_len, column| max(max_len, column.len()));

        let mut result = String::new();

        if r > 100 {
            let lc1 = ((r as f64).log10() as usize) + 5;
            result.push_str(&tab("", lc1));

            let mut space_map: IndexMap<String, usize> = IndexMap::new();
            for k in self.data.keys() {
                let v = &self[&k];
                let mut space = 0usize;
                for elem in v.clone().into_iter().take(5) {
                    space = max(
                        space,
                        min(format!("{:.4}", elem).len(), elem.to_string().len()),
                    );
                }
                if v.len() >= r - 5 {
                    for elem in v.into_iter().skip(r - 5) {
                        space = max(
                            space,
                            min(format!("{:.4}", elem).len(), elem.to_string().len()),
                        );
                    }
                }
                space = max(space + 1, 5);
                if k.len() >= space {
                    space = k.len() + 1;
                }
                result.push_str(&tab(k, space));
                space_map.insert(k.to_string(), space);
            }
            result.push('\n');

            for i in 0..5 {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for k in self.data.keys() {
                    let v = &self[&k];
                    let space = space_map[k];
                    if i < v.len() {
                        let elem = v[i];
                        let st1 = format!("{:.4}", elem);
                        let st2 = elem.to_string();
                        let mut st = st2.clone();

                        if st1.len() < st2.len() {
                            st = st1;
                        }

                        result.push_str(&tab(&st, space));
                    } else {
                        result.push_str(&tab("", space));
                    }
                }
                result.push('\n');
            }
            result.push_str(&tab("...", lc1));
            for k in self.data.keys() {
                let space = space_map[k];
                result.push_str(&tab("...", space));
            }
            result.push('\n');
            for i in r - 5..r {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for k in self.data.keys() {
                    let v = &self[&k];
                    let space = space_map[k];
                    if i < v.len() {
                        let elem = v[i];
                        let st1 = format!("{:.4}", elem);
                        let st2 = elem.to_string();
                        let mut st = st2.clone();

                        if st1.len() < st2.len() {
                            st = st1;
                        }

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

        let mut space_map: IndexMap<String, usize> = IndexMap::new();

        for k in self.data.keys() {
            let value = &self[&k];
            let mut space = 0usize;
            for elem in value {
                space = max(
                    space,
                    min(format!("{:.4}", elem).len(), elem.to_string().len()),
                );
            }
            space = max(space + 1, 5);
            if k.len() >= space {
                space = k.len() + 1;
            }
            result.push_str(&tab(k, space));
            space_map.insert(k.to_string(), space);
        }
        result.push('\n');

        for i in 0..r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for k in self.data.keys() {
                let v = &self[&k];
                let space = space_map[k];
                if i < v.len() {
                    let elem = v[i];
                    let st1 = format!("{:.4}", elem);
                    let st2 = elem.to_string();
                    let mut st = st2.clone();

                    if st1.len() < st2.len() {
                        st = st1;
                    }

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

    /// Pandas-like head (Only print)
    ///
    /// # Description
    /// Print first n rows
    pub fn head_print(&self, n: usize) {
        let r: usize = self.data.values().fold(0, |max_len, column| max(max_len, column.len()));
        let r = n.min(r);
        let mut result = String::new();
        
        result.push_str(&tab("", 5));

        let mut space_map: IndexMap<String, usize> = IndexMap::new();

        for k in self.data.keys() {
            let value = &self[&k];
            let mut space = 0usize;
            for elem in value {
                space = max(
                    space,
                    min(format!("{:.4}", elem).len(), elem.to_string().len()),
                );
            }
            space = max(space + 1, 5);
            if k.len() >= space {
                space = k.len() + 1;
            }
            result.push_str(&tab(k, space));
            space_map.insert(k.to_string(), space);
        }
        result.push('\n');

        for i in 0..r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for k in self.data.keys() {
                let v = &self[&k];
                let space = space_map[k];
                if i < v.len() {
                    let elem = v[i];
                    let st1 = format!("{:.4}", elem);
                    let st2 = elem.to_string();
                    let mut st = st2.clone();

                    if st1.len() < st2.len() {
                        st = st1;
                    }

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
        println!("{}", result);
    }

    /// Pandas-like tail (Only print)
    ///
    /// # Description
    /// Print last n rows
    pub fn tail_print(&self, n: usize) {
        let r: usize = self.data.values().fold(0, |max_len, column| max(max_len, column.len()));
        let mut result = String::new();

        let lc1 = ((r as f64).log10() as usize) + 5;
        result.push_str(&tab("", lc1));

        let mut space_map: IndexMap<String, usize> = IndexMap::new();

        for k in self.data.keys() {
            let value = &self[&k];
            let mut space = 0usize;
            for elem in value {
                space = max(
                    space,
                    min(format!("{:.4}", elem).len(), elem.to_string().len()),
                );
            }
            space = max(space + 1, 5);
            if k.len() >= space {
                space = k.len() + 1;
            }
            result.push_str(&tab(k, space));
            space_map.insert(k.to_string(), space);
        }
        result.push('\n');

        for i in r-n..r {
            result.push_str(&tab(&format!("r[{}]", i), lc1));
            for k in self.data.keys() {
                let v = &self[&k];
                let space = space_map[k];
                if i < v.len() {
                    let elem = v[i];
                    let st1 = format!("{:.4}", elem);
                    let st2 = elem.to_string();
                    let mut st = st2.clone();

                    if st1.len() < st2.len() {
                        st = st1;
                    }

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
        println!("{}", result);
    }
}

/// To deal with CSV file format
pub trait WithCSV: Sized {
    fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>>;
    fn read_csv(file_path: &str, delimiter: char) -> Result<Self, Box<dyn Error>>;
}

/// CSV with DataFrame (inefficient)
impl WithCSV for DataFrame {
    fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r: usize = self
            .data
            .values()
            .fold(0, |max_len, column| max(max_len, column.len()));
        let c: usize = self.data.len();
        wtr.write_record(
            self.data
                .keys()
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
        )?;

        for i in 0..r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for (j, v) in self.data.values().enumerate() {
                if i < v.len() {
                    record[j] = v[i].to_string();
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
        let mut result = DataFrame::with_header(headers);

        for rec in rdr.deserialize() {
            let record: HashMap<String, String> = rec?;
            for head in record.keys() {
                let value = &record[head];
                if value.len() > 0 {
                    (&mut result[&head]).push(value.parse::<f64>().unwrap());
                }
            }
        }

        Ok(result)
    }
}

/// To deal with NetCDF file format
pub trait WithNetCDF: Sized {
    fn write_nc(&self, file_path: &str) -> Result<(), Box<dyn Error>>;
    fn read_nc(file_path: &str) -> Result<Self, Box<dyn Error>>;
    fn read_nc_by_header(file_path: &str, header: Vec<&str>) -> Result<Self, Box<dyn Error>>;
}

/// NetCDF with DataFrame (efficient)
impl WithNetCDF for DataFrame {
    fn write_nc(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut f = netcdf::create(file_path)?;

        for (i, (k, v)) in self.data.iter().enumerate() {
            let dim_name = format!("{}th col", i);
            let dim = v.len();
            f.add_dimension(&dim_name, dim)?;
            let var = &mut f.add_variable::<f64>(k, &[&dim_name])?;
            var.put_values(v, None, None)?;
        }

        Ok(())
    }

    fn read_nc(file_path: &str) -> Result<Self, Box<dyn Error>> {
        let f = netcdf::open(file_path)?;
        let mut df = DataFrame::new();
        for v in f.variables() {
            let k = v.name();
            let mut data: Vec<f64> = vec![0.0; v.len()];
            v.values_to(&mut data, None, None)?;
            df.insert(&k, data);
        }
        Ok(df)
    }

    fn read_nc_by_header(file_path: &str, header: Vec<&str>) -> Result<Self, Box<dyn Error>> {
        let f = netcdf::open(file_path)?;
        let mut df = DataFrame::with_header(header.clone());
        for k in header {
            let val = match f.variable(k) {
                Some(v) => v,
                None => panic!("There are no corresponding values"),
            };
            let mut data: Vec<f64> = vec![0.0; val.len()];
            val.values_to(&mut data, None, None)?;
            df[k] = data;
        }
        Ok(df)
    }
}

pub trait WithJSON {
    fn to_json_value(&self) -> JsonValue;
    fn from_json_value(val: JsonValue) -> Self;
}

impl WithJSON for DataFrame {
    fn to_json_value(&self) -> JsonValue {
        let r = self
            .data
            .values()
            .fold(0, |max_len, column| max(max_len, column.len()));
        let mut values = Vec::<JsonValue>::new();
        for i in 0..r {
            let mut row_object = JsonValue::new_object();
            for head in self.headers() {
                row_object
                    .insert(head, self.data[head][i])
                    .expect("Can't insert row object");
            }
            values.push(row_object);
        }

        JsonValue::Array(values)
    }

    fn from_json_value(val: JsonValue) -> Self {
        unimplemented!()
    }
}
