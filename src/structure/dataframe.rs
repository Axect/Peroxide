extern crate indexmap;
extern crate csv;
extern crate netcdf;

use indexmap::IndexMap;
use self::csv::{ ReaderBuilder, StringRecord, WriterBuilder };
use std::error::Error;
use std::ops::{ Index, IndexMut };
use std::{ fmt, hash::Hash, fmt::Debug };
use std::collections::HashMap;
use structure::matrix::{ Matrix, Shape::*, matrix };
use std::cmp::{ max, min };
use util::useful::tab;

#[derive(Debug, Clone)]
pub struct DataFrame {
    pub data: IndexMap<String, Vec<f64>>,
}

impl fmt::Display for DataFrame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

impl Index<&str> for DataFrame {
    type Output = Vec<f64>;

    fn index(&self, index: &str) -> &Self::Output {
        self.get(index)
    }
}

impl IndexMut<&str> for DataFrame {
    fn index_mut(&mut self, index: &str) -> &mut Self::Output {
        self.data.get_mut(index).unwrap()
    }
}

#[allow(unused_parens)]
impl DataFrame {
    pub fn new() -> Self {
        DataFrame {
            data: IndexMap::new()
        }
    }

    pub fn with_header(header: Vec<&str>) -> Self {
        let l = header.len();
        Self::from_matrix(
            header,
            Matrix {
                data: vec![],
                row: 0,
                col: l,
                shape: Col,
            }
        )
    }

    pub fn insert(&mut self, key: &str, value: Vec<f64>) {
        self.data.insert(key.to_owned(), value);
    }

    pub fn insert_row(&mut self, value: Vec<f64>) {
        assert_eq!(self.data.len(), value.len());
        for (v, val) in self.data.values_mut().zip(value) {
            v.push(val);
        }
    }

    pub fn get(&self, head: &str) -> &Vec<f64> {
        &self.data.get(head).unwrap()
    }

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

    pub fn from_matrix(header: Vec<&str>, mat: Matrix) -> Self {
        let mut df: DataFrame = DataFrame::new();
        for i in 0 .. mat.col {
            df.insert(header[i], mat.col(i));
        }
        df
    }

    pub fn spread(&self) -> String {
        let r: usize = self.data.values().fold(0, |val, v2| max(val, v2.len()));

        let mut result = String::new();

        if r > 100 {
            let lc1 = ((r as f64).log10() as usize) + 5;
            result.push_str(&tab("", lc1));
            let mut space: usize = {
                let mut l = 0usize;
                for v in self.data.values() {
                    for elem in v.clone().into_iter().take(5) {
                        let l2 = min(format!("{:.4}", elem).len(), elem.to_string().len());
                        l = max(l, l2);
                    }
                    if v.len() < r-5 {
                        continue
                    } else {
                        for elem in v.into_iter().skip(r-5) {
                            let l2 = min(format!("{:.4}", elem).len(), elem.to_string().len());
                            l = max(l, l2);
                        }

                    }
                }
                l + 1
            };

            if space < 5 {
                space = 5;
            }
            
            for k in self.data.keys() {
                result.push_str(&tab(k, space));
            }
            result.push('\n');

            for i in 0 .. 5 {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for v in self.data.values() {
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
            for _k in self.data.keys() {
                result.push_str(&tab("...", space));
            }
            result.push('\n');
            for i in r-5 .. r {
                result.push_str(&tab(&format!("r[{}]", i), lc1));
                for v in self.data.values() {
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
                if i == r-1 {
                    break;
                }
                result.push('\n');
            }
            return result;
        }

        result.push_str(&tab("", 5));
        let mut space: usize = {
            let mut l = 0usize;
            for v in self.data.values() {
                for elem in v.into_iter() {
                    let l2 = min(format!("{:.4}", elem).len(), elem.to_string().len());
                    l = max(l, l2);
                }
            }
            l + 1
        };

        if space < 5 {
            space = 5;
        }

        for k in self.data.keys() {
            result.push_str(&tab(k, space));
        }
        result.push('\n');

        for i in 0 .. r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for v in self.data.values() {
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
    
}

pub trait WithCSV: Sized {
    fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>>;
    fn read_csv(file_path: &str, delimiter: char) -> Result<Self, Box<dyn Error>>;
}

impl WithCSV for DataFrame {
    fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r: usize = self.data.values().fold(0, |val, v2| max(val, v2.len()));
        let c: usize = self.data.len();
        wtr.write_record(self.data.keys().map(|x| x.to_string()).collect::<Vec<String>>())?;

        for i in 0 .. r {
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
        let l = headers.len();
        let mut result = DataFrame::from_matrix(headers, Matrix {
            data: vec![], 
            row: 0, 
            col: l, 
            shape: Col
        });

        for rec in rdr.deserialize() {
            let record: HashMap<String, String> = rec?;
            println!("{:?}", record);
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

pub trait WithNetCDF: Sized { 
    fn write_cdf(&self, file_path: &str) -> Result<(), Box<dyn Error>>;         
    fn read_cdf(file_path: &str, header: Vec<&str>) -> Result<Self, Box<dyn Error>>;
}                                                                               
                                                                                
impl WithNetCDF for DataFrame {         
    fn write_cdf(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut f = netcdf::create(file_path)?;
        
        for (i, (k, v)) in self.data.iter().enumerate() {
            let dim_name = format!("{}th col", i);
            let dim = v.len();
            f.root.add_dimension(&dim_name, dim as u64)?;
            f.root.add_variable(
                k,
                &vec![dim_name],
                v
            )?;
        }

        Ok(())
    }
    fn read_cdf(file_path: &str, header: Vec<&str>) -> Result<Self, Box<dyn Error>> {
        let f = netcdf::open(file_path)?;
        
        unimplemented!()
    }
}
