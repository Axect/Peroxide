extern crate indexmap;
use indexmap::IndexMap;
use std::ops::{ Index, IndexMut };
use std::{ fmt, hash::Hash, fmt::Debug };
use structure::matrix::{ Matrix, Shape::*, matrix };
use std::cmp::{ max, min };
use util::useful::tab;

#[derive(Debug, Clone)]
pub struct DataFrame<T> where T: Hash + Eq + Clone {
    pub data: IndexMap<T, Vec<f64>>,
}

impl<T> fmt::Display for DataFrame<T> where T: Hash + Eq + Clone + Debug {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

impl<T> Index<T> for DataFrame<T> where T: Hash + Clone + Eq + Debug {
    type Output = Vec<f64>;

    fn index(&self, index: T) -> &Self::Output {
        self.get(index)
    }
}

impl<T> IndexMut<T> for DataFrame<T> where T: Hash + Clone + Eq + Debug {
    fn index_mut(&mut self, index: T) -> &mut Self::Output {
        self.data.get_mut(&index).unwrap()
    }
}

#[allow(unused_parens)]
impl<T> DataFrame<T> where T: Hash + Eq + Clone + Debug {
    pub fn new() -> Self {
        DataFrame {
            data: IndexMap::new()
        }
    }

    pub fn insert(&mut self, key: T, value: Vec<f64>) {
        self.data.insert(key, value);
    }

    pub fn get(&self, head: T) -> &Vec<f64> {
        &self.data.get(&head).unwrap()
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

    pub fn from_matrix(header: Vec<T>, mat: Matrix) -> Self {
        let mut df: DataFrame<T> = DataFrame::new();
        for i in 0 .. mat.col {
            df.insert(header[i].clone(), mat.col(i));
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
                    let v_len = v.len();
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
                result.push_str(&tab(&format!("{:?}", k), space));
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
            result.push_str(&tab(&format!("{:?}", k), space));
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

