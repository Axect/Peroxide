extern crate csv;

#[allow(unused_imports)]
use structure::vector::*;
use std::convert;
use std::fmt;
use std::ops::{Add, Neg, Sub, Mul, Rem, Index, IndexMut};
pub use self::Shape::{Row, Col};
pub use self::Norm::*;
use std::cmp::{max, min};
pub use std::error::Error;
use self::csv::{WriterBuilder, ReaderBuilder, StringRecord};
use util::useful::*;

pub type Perms = Vec<(usize, usize)>;

/// To select matrices' binding.
/// 
/// Row - Row binding
/// Col - Column binding
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4], 2, 2, Row);
/// let b = matrix(vec![1,2,3,4], 2, 2, Col);
/// println!("{}", a); // Similar to [[1,2],[3,4]]
/// println!("{}", b); // Similar to [[1,3],[2,4]]
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Shape {
    Row,
    Col,
}

/// Print for Shape
impl fmt::Display for Shape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let to_display = match self {
            Row => "Row",
            Col => "Col",
        };
        write!(f, "{}", to_display)
    }
}

/// R-like matrix structure
///
/// # Examples
///
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = Matrix {
///     data: vec![1f64,2f64,3f64,4f64],
///     row: 2,
///     col: 2,
///     shape: Row,
/// }; // [[1,2],[3,4]]
/// ```
#[derive(Debug, Clone)]
pub struct Matrix {
    pub data: Vec<f64>,
    pub row: usize,
    pub col: usize,
    pub shape: Shape,
}

// =============================================================================
// Various matrix constructor
// =============================================================================

/// R-like matrix constructor
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(c!(1,2,3,4), 2, 2, Row);
/// ```
pub fn matrix<T>(v: Vec<T>, r: usize, c:usize, shape: Shape) -> Matrix where T: convert::Into<f64> {
    Matrix {
        data: v.into_iter().map(|t| t.into()).collect::<Vec<f64>>(),
        row: r,
        col: c,
        shape
    }
}

/// R-like matrix constructor (Explicit ver.)
pub fn r_matrix<T>(v: Vec<T>, r: usize, c: usize, shape: Shape) -> Matrix where T: convert::Into<f64> {
    matrix(v, r, c, shape)
}

/// Python-like matrix constructor
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = py_matrix(vec![c!(1,2), c!(3,4)]);
/// let b = matrix(c!(1,2,3,4), 2, 2, Row);
/// assert_eq!(a, b);
/// ```
pub fn py_matrix<T>(v: Vec<Vec<T>>) -> Matrix where T: convert::Into<f64> + Copy {
    let r = v.len();
    let c = v[0].len();
    let mut data = vec![0f64; r*c];
    for i in 0 .. r {
        for j in 0 .. c {
            let idx = i * c + j;
            data[idx] = v[i][j].into();
        }
    }
    matrix(data, r, c, Row)
}

/// Matlab-like matrix constructor
/// 
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
/// 
/// let a = ml_matrix("1 2; 3 4");
/// let b = matrix(c!(1,2,3,4), 2, 2, Row);
/// assert_eq!(a, b);
/// ```
pub fn ml_matrix(s: &str) -> Matrix where {
    let str_rows: Vec<&str> = s.split(';').collect();
    let r = str_rows.len();
    let str_data = str_rows.into_iter()
        .map(|x| x.trim().split(' ').collect::<Vec<&str>>())
        .collect::<Vec<Vec<&str>>>();
    let c = str_data[0].len();
    let data = str_data.into_iter()
        .flat_map(|t| t.into_iter().map(|x| x.parse::<f64>().unwrap()).collect::<Vec<f64>>())
        .collect::<Vec<f64>>();
    matrix(data, r, c, Row)
}

/// Pretty Print
impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

/// PartialEq implements
impl PartialEq for Matrix {
    fn eq(&self, other: &Matrix) -> bool {
        if self.shape == other.shape {
            self.data.clone()
                .into_iter()
                .zip(other.data.clone())
                .all(|(x, y)| nearly_eq(x,y)) && self.row == other.row
        } else {
            self.eq(&other.change_shape())
        }
    }
}

/// Main matrix structure
#[allow(dead_code)]
impl Matrix {
    /// Change Bindings
    ///
    /// `Row` -> `Col` or `Col` -> `Row`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(vec![1,2,3,4],2,2,Row);
    /// assert_eq!(a.shape, Row);
    /// let b = a.change_shape();
    /// assert_eq!(b.shape, Col);
    /// ```
    pub fn change_shape(&self) -> Matrix {
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

    /// Spread data(1D vector) to 2D formatted String
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(vec![1,2,3,4],2,2,Row);
    /// println!("{}", a.spread()); // same as println!("{}", a);
    /// // Result:
    /// //       c[0] c[1]
    /// // r[0]     1    3
    /// // r[1]     2    4
    /// ```
    pub fn spread(&self) -> String {
        assert_eq!(self.row * self.col, self.data.len());
        let r = self.row;
        let c = self.col;
        let mut key_row = 20usize;
        let mut key_col = 20usize;

        if r > 100 || c > 100 || (r > 20 && c > 20) {
            let part =
                if r <= 10 {
                    key_row = r;
                    key_col = 100;
                    self.take(100, Col)
                } else if c <= 10 {
                    key_row = 100;
                    key_col = c;
                    self.take(100, Row)
                } else {
                    self.take(20, Row).take(20, Col)
                };
            return format!(
                "Result is too large to print - {}x{}\nonly print {}x{} parts:\n{}",
                self.row.to_string(),
                self.col.to_string(),
                key_row.to_string(),
                key_col.to_string(),
                part.spread()
            );
        }

        // Find maximum length of data
        let sample = self.data.clone();
        let mut space: usize = sample.into_iter()
            .map(|x|
                min(
                    format!("{:.4}", x).len(),
                    x.to_string().len(),
                ) // Choose minimum of approx vs normal
            )
            .fold(0, |x, y| max(x,y)) + 1;

        if space < 5 {
            space = 5;
        }

        let mut result = String::new();

        result.push_str(&tab("", 5));
        for i in 0 .. c {
            result.push_str(&tab(&format!("c[{}]", i), space)); // Header
        }
        result.push('\n');

        for i in 0 .. r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for j in 0 .. c {
                let st1 = format!("{:.4}",self[(i, j)]); // Round at fourth position
                let st2 = self[(i,j)].to_string();       // Normal string
                let mut st = st2.clone();

                // Select more small thing
                if st1.len() < st2.len() {
                    st = st1;
                }

                result.push_str(&tab(&st, space));
            }
            if i == (r-1) {
                break;
            }
            result.push('\n');
        }

        return result;
    }

    /// Extract Column
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(c!(1,2,3,4), 2, 2, Row);
    /// assert_eq!(a.col(0), c!(1,3));
    /// ```
    pub fn col(&self, index: usize) -> Vector {
        assert!(index < self.col);
        let mut container: Vec<f64> = Vec::new();
        match self.shape {
            Row => {
                let l: usize = self.row * self.col;
                for i in 0 .. l {
                    if i % self.col == index {
                        container.push(self.data[i]);
                    }
                }
            },
            Col => {
                let s: usize = self.row * index;
                container = self.data.clone().into_iter()
                    .skip(s)
                    .take(self.row).collect::<Vec<f64>>();
            }
        }
        container
    }

    /// Extract Row
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(c!(1,2,3,4), 2, 2, Row);
    /// assert_eq!(a.row(0), c!(1,2));
    /// ```
    pub fn row(&self, index: usize) -> Vector {
        assert!(index < self.row);
        let mut container: Vec<f64> = Vec::new();
        match self.shape {
            Row => {
                let s: usize = self.col * index;
                container = self.data.clone().into_iter()
                    .skip(s)
                    .take(self.col).collect::<Vec<f64>>();
            },
            Col => {
                let l: usize = self.row * self.col;
                for i in 0 .. l {
                    if i % self.row == index {
                        container.push(self.data[i]);
                    }
                }
            }
        }
        container
    }

    /// Extract diagonal components
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.diag(), c!(1,4));
    /// ```
    pub fn diag(&self) -> Vector {
        let mut container: Vector = Vec::new();
        let r = self.row;
        let c = self.col;
        assert_eq!(r, c);
        for i in 0 .. r {
            container.push(self.data[i * (c + 1)]);
        }
        container
    }

    /// Swap row or col
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.swap(0,1,Row), matrix(c!(3,4,1,2),2,2,Row));
    /// assert_eq!(a.swap(0,1,Col), matrix(c!(2,4,1,3),2,2,Col));
    /// ```
    pub fn swap(&self, idx1: usize, idx2: usize, shape: Shape) -> Matrix {
        match shape {
            Row => {
                let mut v: Vector = Vec::new();
                for k in 0 .. self.row {
                    if k == idx1 {
                        v.extend(&self.row(idx2));
                    } else if k == idx2 {
                        v.extend(&self.row(idx1));
                    } else {
                        v.extend(&self.row(k));
                    }
                }
                matrix(v, self.row, self.col, Row)
            },
            Col => {
                let mut v: Vector = Vec::new();
                for k in 0 .. self.col {
                    if k == idx1 {
                        v.extend(&self.col(idx2));
                    } else if k == idx2 {
                        v.extend(&self.col(idx1));
                    } else {
                        v.extend(&self.col(k));
                    }
                }
                matrix(v, self.row, self.col, Col)
            }
        }
    }

    /// Write to CSV
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// a.write("test.csv");
    /// ```
    pub fn write(&self, file_path: &str) -> Result<(), Box<Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r = self.row;
        let c = self.col;
        for i in 0 .. r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for j in 0 .. c {
                record[j] = self[(i, j)].to_string();
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    /// Write to CSV (with round option)
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// a.write_round("test.csv", 0);
    /// ```
    pub fn write_round(&self, file_path: &str, round: usize) -> Result<(), Box<Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r = self.row;
        let c = self.col;
        for i in 0 .. r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for j in 0 .. c {
                record[j] = format!("{:.*}", round, self[(i, j)]);
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    pub fn write_with_header(&self, file_path: &str, header: Vec<&str>) -> Result<(), Box<Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r = self.row;
        let c = self.col;
        assert_eq!(c, header.len());
        wtr.write_record(header)?;
        for i in 0 .. r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for j in 0 .. c {
                record[j] = self[(i, j)].to_string();
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    pub fn write_with_header_round(&self, file_path: &str, header: Vec<&str>, round: usize) -> Result<(), Box<Error>> {
        let mut wtr = WriterBuilder::new().from_path(file_path)?;
        let r = self.row;
        let c = self.col;
        assert_eq!(c, header.len());
        wtr.write_record(header)?;
        for i in 0 .. r {
            let mut record: Vec<String> = vec!["".to_string(); c];
            for j in 0 .. c {
                record[j] = format!("{:.*}", round, self[(i, j)]);
            }
            wtr.write_record(record)?;
        }
        wtr.flush()?;
        Ok(())
    }

    /// Read from CSV
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    /// use std::process;
    ///
    /// let a = matrix(c!(1,2,3,3,2,1), 3, 2, Col);
    /// a.write_round("test.csv", 0);
    ///
    /// let b = Matrix::read("test.csv", false, ','); // header = false, delimiter = ','
    /// match b {
    ///     Ok(mat) => println!("{}", mat),
    ///     Err(err) => {
    ///         println!("{}", err);
    ///         process::exit(1);
    ///     }
    /// }
    /// ```
    pub fn read(file_path: &str, header: bool, delimiter: char) -> Result<Matrix, Box<Error>> {
        let mut rdr = ReaderBuilder::new().has_headers(header)
            .delimiter(delimiter as u8)
            .from_path(file_path)?;

        let records = rdr.records().collect::<Result<Vec<StringRecord>, csv::Error>>()?;
        let result = records;

        let l_row = result.len();
        let l_col = result[0].len();

        let mut m = matrix(vec![0f64; l_row * l_col], l_row, l_col, Row);

        for i in 0 .. l_row {
            for j in 0 .. l_col {
                m[(i, j)] = result[i][j].parse::<f64>().unwrap();
            }
        }

        Ok(m)
    }

    /// Substitute Col
    pub fn subs_col(&mut self, idx: usize, v: Vec<f64>) {
        for i in 0 .. self.row {
           self[(i, idx)] = v[i]; 
        }
    }

    /// Substitute Row
    pub fn subs_row(&mut self, idx: usize, v: Vec<f64>) {
        for j in 0 .. self.col {
            self[(idx, j)] = v[j];
        }
    }

    /// From index operations
    pub fn from_index<F, G>(f: F, size: (usize, usize)) -> Matrix
    where F: Fn(usize, usize) -> G + Copy,
          G: Into<f64>
    {
        let row = size.0;
        let col = size.1;

        let mut mat = matrix(vec![0f64; row*col], row, col, Row);

        for i in 0 .. row {
            for j in 0 .. col {
                mat[(i, j)] = f(i, j).into();
            }
        }
        mat
    }
}


// =============================================================================
// Common Properties of Matrix & Vector
// =============================================================================

/// Common trait for Matrix & Vector
pub trait LinearOps {
    fn to_matrix(&self) -> Matrix;
    fn transpose(&self) -> Matrix;
    fn t(&self) -> Matrix;
}

impl LinearOps for Matrix {
    /// Just clone
    fn to_matrix(&self) -> Matrix {
        self.clone()
    }

    /// Transpose
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(vec![1,2,3,4], 2, 2, Row);
    /// println!("{}", a); // [[1,3],[2,4]]
    /// ```
    fn transpose(&self) -> Matrix {
        match self.shape {
            Row => matrix(self.data.clone(), self.col, self.row, Col),
            Col => matrix(self.data.clone(), self.col, self.row, Row)
        }
    }

    /// R-like transpose function
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.transpose(), a.t());
    /// ```
    fn t(&self) -> Matrix {
        self.transpose()
    }
}

impl LinearOps for Vector {
    /// Vector to Column Matrix
    fn to_matrix(&self) -> Matrix {
        let l = self.len();
        matrix(self.clone(), l, 1, Col)
    }

    /// Vector to Row Matrix
    fn transpose(&self) -> Matrix {
        let l = self.len();
        matrix(self.clone(), 1, l, Row)
    }

    /// R-like Syntax
    fn t(&self) -> Matrix {
        self.transpose()
    }
}

// =============================================================================
// Standard Operation for Matrix (ADD)
// =============================================================================

/// Element-wise addition of Matrix
///
/// # Caution
/// > You should remember ownership.
/// > If you use Matrix `a,b` then you can't use them after.
impl Add<Matrix> for Matrix {
    type Output = Matrix;

    fn add(self, other: Matrix) -> Matrix {
        assert_eq!(&self.row, &other.row);
        assert_eq!(&self.col, &other.col);
        if self.shape == other.shape {
            matrix(
                self.data.clone().into_iter().zip(&other.data).map(|(x,y)| x + y).collect::<Vec<f64>>(),
                self.row,
                self.col,
                self.shape,
            )
        } else {
            self.change_shape().add(other)
        }
    }
}

/// Element-wise addition between Matrix & f64
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// assert_eq!(a + 1, matrix!(2;5;1, 2, 2, Row));
/// ```
impl<T> Add<T> for Matrix where T: convert::Into<f64> + Copy {
    type Output = Matrix;
    fn add(self, other: T) -> Matrix {
        self.fmap(|x| x + other.into())
    }
}

/// Element-wise addition between f64 & matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// assert_eq!(1f64 + a, matrix!(2;5;1, 2, 2, Row));
/// ```
impl Add<Matrix> for f64 {
    type Output = Matrix;

    fn add(self, other: Matrix) -> Matrix {
        other.add(self)
    }
}

/// Element-wise addition between f64 & matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// assert_eq!(1 + a, matrix!(2;5;1, 2, 2, Row));
/// ```
impl Add<Matrix> for i32 {
    type Output = Matrix;

    fn add(self, other: Matrix) -> Matrix {
        other.add(self)
    }
}

/// Element-wise addition between f64 & matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// assert_eq!(1 as usize + a, matrix!(2;5;1, 2, 2, Row));
/// ```
impl Add<Matrix> for usize {
    type Output = Matrix;

    fn add(self, other: Matrix) -> Matrix {
        other.add(self as f64)
    }
}


// =============================================================================
// Standard Operation for Matrix (Neg)
// =============================================================================
/// Negation of Matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4],2,2,Row);
/// println!("{}", -a); // [[-1,-2],[-3,-4]]
/// ```
impl Neg for Matrix {
    type Output = Matrix;
    
    fn neg(self) -> Matrix {
        matrix(
            self.data.into_iter().map(|x:f64| -x).collect::<Vec<f64>>(),
            self.row,
            self.col,
            self.shape,
        )
    }
}

// =============================================================================
// Standard Operation for Matrix (Sub)
// =============================================================================
/// Subtraction between Matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4],2,2,Row);
/// let b = matrix(vec![1,2,3,4],2,2,Col);
/// println!("{}", a - b); // [[0, -1], [1, 0]]
/// ```
impl Sub<Matrix> for Matrix {
    type Output = Matrix;

    fn sub(self, other: Matrix) -> Matrix {
        assert_eq!(&self.row, &other.row);
        assert_eq!(&self.col, &other.col);
        if self.shape == other.shape {
            matrix(
                self.data.clone().into_iter().zip(&other.data).map(|(x,y)| x - y).collect::<Vec<f64>>(),
                self.row,
                self.col,
                self.shape,
            )
        } else {
            self.change_shape().sub(other)
        }
    }
}

impl<T> Sub<T> for Matrix where T: convert::Into<f64> + Copy {
    type Output = Matrix;

    fn sub(self, other: T) -> Matrix {
        self.fmap(|x| x - other.into())
    }
}

/// Subtraction Matrix with f64
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4],2,2,Row);
/// assert_eq!(a - 1f64, matrix!(0;3;1, 2, 2, Row));
/// ```
impl Sub<Matrix> for f64 {
    type Output = Matrix;

    fn sub(self, other: Matrix) -> Matrix {
        - other.sub(self)
    }
}

impl Sub<Matrix> for i32 {
    type Output = Matrix;

    fn sub(self, other: Matrix) -> Matrix {
        - other.sub(self)
    }
}

impl Sub<Matrix> for usize {
    type Output = Matrix;

    fn sub(self, other: Matrix) -> Matrix {
        - other.sub(self as f64)
    }
}

/// Element-wise matrix multiplication
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4], 2, 2, Row);
/// let b = matrix(vec![1,2,3,4], 2, 2, Col);
/// println!("{}", a * b); // [[1,6],[6,16]]
/// ```
impl Mul<Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        assert_eq!(self.row, other.row);
        assert_eq!(self.col, other.col);
        self.zip_with(|x,y| x * y, &other)
    }
}

impl<T> Mul<T> for Matrix where T: convert::Into<f64> + Copy {
    type Output = Matrix;

    fn mul(self, other: T) -> Matrix {
        self.fmap(|x| x * other.into())
    }
}

impl Mul<Matrix> for f64 {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        other.mul(self)
    }
}

impl Mul<Matrix> for i64 {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        other.mul(self as f64)
    }
}

impl Mul<Matrix> for i32 {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        other.mul(self)
    }
}

impl Mul<Matrix> for usize {
    type Output = Matrix;

    fn mul(self, other: Matrix) -> Matrix {
        other.mul(self as f64)
    }
}

/// Matrix Multiplication
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// let b = matrix!(1;4;1, 2, 2, Col);
/// assert_eq!(a % b, matrix(c!(5, 11, 11, 25), 2, 2, Row));
///
/// let m = matrix!(1;4;1, 2, 2, Row);
/// let v = c!(1,2);
/// assert_eq!(m % v, matrix(c!(5,11),2,1,Col));
/// ```
impl<T> Rem<T> for Matrix where T: LinearOps {
    type Output = Matrix;

    fn rem(self, other: T) -> Matrix {
        let r_self = self.row;
        let c_self = self.col;
        let new_other = other.to_matrix();
        let r_other = new_other.row;
        let c_other = new_other.col;

        assert_eq!(c_self, r_other);

        let r_new = r_self;
        let c_new = c_other;

        let mut result = matrix(vec![0f64; r_new * c_new], r_new, c_new, self.shape);

        for i in 0 .. r_new {
            for j in 0 .. c_new {
                let mut s = 0f64;
                for k in 0 .. c_self {
                    s += self[(i, k)] * new_other[(k, j)];
                }
                result[(i, j)] = s;
            }
        }
        result
    }
}

/// Matrix multiplication for Vector vs Matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;4;1, 2, 2, Row);
/// let v = c!(1,2);
/// assert_eq!(v % a, matrix(c!(7,10),1,2,Row));
/// ```
impl Rem<Matrix> for Vector {
    type Output = Matrix;

    fn rem(self, other: Matrix) -> Matrix {
        assert_eq!(self.len(), other.row);
        let l = self.len();
        matrix(self, 1, l, Col).rem(other)
    }
}

/// Index for Matrix
///
/// `(usize, usize) -> f64`
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(vec![1,2,3,4],2,2,Row);
/// assert_eq!(a[(0,1)], 2f64);
/// ```
impl Index<(usize, usize)> for Matrix {
    type Output = f64;

    fn index(&self, pair: (usize, usize)) -> &f64 {
        let i = pair.0;
        let j = pair.1;
        match self.shape {
            Row => &self.data[i * self.col + j],
            Col => &self.data[i + j * self.row]
        }
    }
}

/// IndexMut for Matrix (Assign)
///
/// `(usize, usize) -> f64`
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let mut a = matrix!(1;4;1, 2, 2, Row);
/// a[(1,1)] = 10.0;
/// assert_eq!(a, matrix(c!(1,2,3,10), 2, 2, Row));
/// ```
impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, pair: (usize, usize)) -> &mut f64 {
        let i = pair.0;
        let j = pair.1;
        let r = self.row;
        let c = self.col;
        match self.shape {
            Row => {
                let idx = i * c + j;
                &mut self.data[idx]
            },
            Col => {
                let idx = i + j * r;
                &mut self.data[idx]
            }
        }
    }
}

// =============================================================================
// Functional Programming Tools (Hand-written)
// =============================================================================
pub trait FP {
    fn take(&self, n: usize, shape: Shape) -> Matrix;
    fn skip(&self, n: usize, shape: Shape) -> Matrix;
    fn fmap<F>(&self, f: F) -> Matrix where F: Fn(f64) -> f64;
    fn reduce<F, T>(&self, init: T, f: F) -> f64 
        where F: Fn(f64, f64) -> f64,
              T: convert::Into<f64>;
    fn zip_with<F>(&self, f: F, other: &Matrix) -> Matrix 
        where F: Fn(f64, f64) -> f64;
}

impl FP for Matrix {
    fn take(&self, n: usize, shape: Shape) -> Matrix {
        match shape {
            Row => {
                let mut temp_data: Vec<f64> = Vec::new();
                for i in 0..n {
                    if i >= self.row {
                        return matrix(temp_data, self.row, self.col, Row);
                    }
                    temp_data.extend(self.row(i));
                }
                matrix(temp_data, n, self.col, Row)
            },
            Col => {
                let mut temp_data: Vec<f64> = Vec::new();
                for i in 0..n {
                    if i >= self.col {
                        return matrix(temp_data, self.row, self.col, Col);
                    }
                    temp_data.extend(self.col(i));
                }
                matrix(temp_data, self.row, n, Col)
            }
        }
    }

    fn skip(&self, n: usize, shape: Shape) -> Matrix {
        match shape {
            Row => {
                assert!(n < self.row, "Skip range is larger than row of matrix");

                let mut temp_data: Vec<f64> = Vec::new();
                for i in n .. self.row {
                    temp_data.extend(self.row(i));
                }
                matrix(temp_data, self.row - n, self.col, Row)
            },
            Col => {
                assert!(n < self.col, "Skip range is larger than col of matrix");

                let mut temp_data: Vec<f64> = Vec::new();
                for i in n .. self.col {
                    temp_data.extend(self.col(i));
                }
                matrix(temp_data, self.row, self.col - n, Col)
            }
        }
    }

    fn fmap<F>(&self, f: F) -> Matrix where F: Fn(f64) -> f64 {
        let result = self.data.clone().into_iter().map(|x| f(x)).collect::<Vec<f64>>();
        matrix(
            result,
            self.row,
            self.col,
            self.shape,
        )
    }

    fn reduce<F, T>(&self, init: T, f: F) -> f64 
        where F: Fn(f64, f64) -> f64,
              T: convert::Into<f64> {
        self.data.clone().into_iter().fold(
            init.into(),
            |x,y| f(x,y),
        )
    }

    fn zip_with<F>(&self, f: F, other: &Matrix) -> Matrix 
        where F: Fn(f64, f64) -> f64 {
        assert_eq!(self.data.len(), other.data.len());
        let mut a = other.clone();
        if self.shape != other.shape {
            a = a.change_shape();
        }
        let result = self.data.clone().into_iter()
            .zip(a.data)
            .map(|(x,y)| f(x,y))
            .collect::<Vec<f64>>();
        matrix(
            result,
            self.row,
            self.col,
            self.shape,
        )
    }
}

// =============================================================================
// Linear Algebra
// =============================================================================

/// Norm Enum
#[derive(Debug, Copy, Clone)]
pub enum Norm {
    Frobenius,
    PQ(usize, usize),
    One,
    Infinity
}

/// Linear algebra trait
pub trait LinearAlgebra {
    fn norm(&self, norm: Norm) -> f64;
    fn lu(&self) -> Option<PQLU>;
    fn det(&self) -> f64;
    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix);
    fn inv(&self) -> Option<Matrix>;
    fn pseudo_inv(&self) -> Option<Matrix>;
}

pub fn diag(n: usize) -> Matrix {
    let mut v: Vec<f64> = vec![0f64; n*n];
    for i in 0 .. n {
        let idx = i * (n + 1);
        v[idx] = 1f64;
    }
    matrix(v, n, n, Row)
}

#[derive(Debug, Clone)]
pub struct PQLU {
    pub p: Perms,
    pub q: Perms,
    pub l: Matrix,
    pub u: Matrix,
}

impl LinearAlgebra for Matrix {
    /// Matrix norm
    ///
    /// # Kinds
    /// * `Frobenius` : Frobenius norm
    /// * `PQ(usize, usize)` : L_pq norm
    /// * `One` : 1-norm
    /// * `Infinity` : Infinity norm
    fn norm(&self, norm: Norm) -> f64 {
        match norm {
            Frobenius => {
                let mut s = 0f64;
                for i in 0 .. self.data.len() {
                    s += self.data[i].powi(2);
                }
                s.sqrt()
            },
            PQ(p, q) => {
                let mut s = 0f64;
                for j in 0 .. self.col {
                    let mut s_row = 0f64;
                    for i in 0 .. self.row {
                        s_row += self[(i, j)].powi(p as i32);
                    }
                    s += s_row.powf(q as f64 / (p as f64));
                }
                s.powf(1f64 / (q as f64))
            },
            One => {
                let mut m = std::f64::MIN;
                let a = match self.shape {
                    Row => self.change_shape(),
                    Col => self.clone(),
                };
                for c in 0 .. a.col {
                    let s = a.col(c).reduce(0f64, |x, y| x + y);
                    if s > m {
                        m = s;
                    }
                }
                m
            },
            Infinity => {
                let mut m = std::f64::MIN;
                let a = match self.shape {
                    Row => self.clone(),
                    Col => self.change_shape(),
                };
                for r in 0 .. a.row {
                    let s = a.row(r).reduce(0f64, |x, y| x + y);
                    if s > m {
                        m = s;
                    }
                }
                m
            }
        }
    }

    /// LU Decomposition Implements
    ///
    /// # Description
    /// It use complete pivoting LU decomposition.
    /// You can get two permutations, and LU matrices.
    ///
    /// # Caution
    /// It returns `Option<PQLU>` - You should unwrap to obtain real value.
    /// `PQLU` has four field - `p`, `q`, `l`, `u`.
    /// `p`, `q` are permutations.
    /// `l`, `u` are matrices.
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix(vec![1,2,3,4], 2, 2, Row);
    /// let pqlu = a.lu().unwrap();
    /// let (p,q,l,u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);
    /// assert_eq!(p, vec![(0,1)]); // swap 0 & 1 (Row)
    /// assert_eq!(q, vec![(0,1)]); // swap 0 & 1 (Col)
    /// assert_eq!(l, matrix(c!(1,0,0.5,1),2,2,Row));
    /// assert_eq!(u, matrix(c!(4,3,0,-0.5),2,2,Row));
    /// ```
    fn lu(&self) -> Option<PQLU> {
        assert_eq!(self.col, self.row);
        let n = self.row;
        let len: usize = n * n;

        let mut l: Matrix = matrix(vec![0f64; len], n, n, self.shape);
        let mut u: Matrix = matrix(vec![0f64; len], n, n, self.shape);

        // ---------------------------------------
        // Pivoting - Complete
        // ---------------------------------------
        // Permutations
        let mut p: Perms = Vec::new();
        let mut q: Perms = Vec::new();

        let mut container = self.clone();

        for k in 0 .. (n-1) {
            // Initialize maximum & Position
            let mut m = 0f64;
            let mut row_idx: usize = k;
            let mut col_idx: usize = k;
            // Find Maximum value & Position
            for i in k .. n {
                for j in k .. n {
                    let temp = container[(i,j)];
                    if temp.abs() > m.abs() {
                        m = temp;
                        row_idx = i;
                        col_idx = j;
                    }
                }
            }
            if k != row_idx {
                container = container.swap(k, row_idx, Row); // Row perm
                p.push((k, row_idx));
            }
            if k != col_idx {
                container = container.swap(k, col_idx, Col); // Col perm
                q.push((k, col_idx));
            }
        }

        // ---------------------------------------
        // Obtain L & U
        // ---------------------------------------
        let reference = container.clone();

        // Initialize U
        for i in 0 .. n {
            u[(0, i)] = reference[(0, i)];
        }

        // Initialize L
        for i in 0 .. n {
            l[(i, i)] = 1f64;
        }

        for i in 0 .. n {
            for k in i .. n {
                let mut s = 0f64;
                for j in 0 .. i {
                    s += l[(i, j)] * u[(j, k)];
                }
                u[(i, k)] = reference[(i, k)] - s;
                // Check non-zero diagonal
                if nearly_eq(u[(i, i)], 0) {
                    return None;
                }
            }

            for k in (i+1) .. n {
                let mut s = 0f64;
                for j in 0 .. i {
                    s += l[(k, j)] * u[(j, i)];
                }
                l[(k, i)] = (reference[(k, i)] - s) / u[(i, i)];
            }
        }

        Some(
            PQLU {
                p,
                q,
                l,
                u,
            }
        )
    }

    /// Determinant
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.det(), -2f64);
    /// ```
    fn det(&self) -> f64 {
        assert_eq!(self.row, self.col);
        match self.lu() {
            None => 0f64,
            Some(pqlu) => {
                let (p, q, _l, u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);

                // sgn of perms
                let sgn_p = 2.0 * (p.len() % 2) as f64 - 1.0;
                let sgn_q = 2.0 * (q.len() % 2) as f64 - 1.0;

                u.diag().reduce(1f64, |x,y| x * y) * sgn_p * sgn_q
            }
        }
    }

    /// Block Partition
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;16;1, 4, 4, Row);
    /// let (m1, m2, m3, m4) = a.block();
    /// assert_eq!(m1, matrix(c!(1,2,5,6), 2, 2, Row));
    /// assert_eq!(m2, matrix(c!(3,4,7,8), 2, 2, Row));
    /// assert_eq!(m3, matrix(c!(9,10,13,14), 2, 2, Row));
    /// assert_eq!(m4, matrix(c!(11,12,15,16), 2, 2, Row));
    ///
    /// let b = matrix!(1;16;1, 4, 4, Col);
    /// let (m1, m2, m3, m4) = b.block();
    /// assert_eq!(m1, matrix(c!(1,2,5,6), 2, 2, Col));
    /// assert_eq!(m3, matrix(c!(3,4,7,8), 2, 2, Col));
    /// assert_eq!(m2, matrix(c!(9,10,13,14), 2, 2, Col));
    /// assert_eq!(m4, matrix(c!(11,12,15,16), 2, 2, Col));
    /// ```
    fn block(&self) -> (Matrix, Matrix, Matrix, Matrix) {
        let r = self.row;
        let c = self.col;
        let l = min(self.row / 2, self.col / 2);
        let r_l = r - l;
        let c_l = c - l;

        let mut m1 = matrix(vec![0f64; l * l], l, l, self.shape);
        let mut m2 = matrix(vec![0f64; l * c_l], l, c_l, self.shape);
        let mut m3 = matrix(vec![0f64; r_l * l], r_l, l, self.shape);
        let mut m4 = matrix(vec![0f64; r_l * c_l], r_l, c_l, self.shape);

        for idx_row in 0 .. r {
            for idx_col in 0..c {
                match (idx_row, idx_col) {
                    (i, j) if (i < l) && (j < l) => {
                        m1[(i, j)] = self[(i, j)];
                    },
                    (i, j) if (i < l) && (j >= l) => {
                        m2[(i, j - l)] = self[(i, j)];
                    },
                    (i, j) if (i >= l) && (j < l) => {
                        m3[(i - l, j)] = self[(i, j)];
                    },
                    (i, j) if (i >= l) && (j >= l) => {
                        m4[(i - l, j - l)] = self[(i, j)];
                    },
                    _ => (),
                }
            }
        }
        (m1, m2, m3, m4)
    }

    /// Inverse of Matrix
    ///
    /// # Caution
    ///
    /// `inv` function returns `Option<Matrix>`
    /// Thus, you should use pattern matching or `unwrap` to obtain inverse.
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// // Non-singular
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// assert_eq!(a.inv().unwrap(), matrix(c!(-2,1,1.5,-0.5),2,2,Row));
    ///
    /// // Singular
    /// let b = matrix!(1;9;1, 3, 3, Row);
    /// assert_eq!(b.inv(), None);
    /// ```
    fn inv(&self) -> Option<Matrix> {
        match self.lu() {
            None => None,
            Some(pqlu) => {
                let (p, q, l, u) = (pqlu.p, pqlu.q, pqlu.l, pqlu.u);
                let mut m = inv_u(u) % inv_l(l);
                for (idx1, idx2) in q.into_iter().rev() {
                    m = m.swap(idx1, idx2, Row);
                }
                for (idx1, idx2) in p.into_iter().rev() {
                    m = m.swap(idx1, idx2, Col);
                }
                Some(m)
            }
        }
    }

    /// Moore-Penrose Pseudo inverse
    ///
    /// # Description
    /// `(X^T X)^{-1} X`
    ///
    /// # Examples
    /// ```
    /// extern crate peroxide;
    /// use peroxide::*;
    ///
    /// let a = matrix!(1;4;1, 2, 2, Row);
    /// let inv_a = a.inv().unwrap();
    /// let pse_a = a.pseudo_inv().unwrap();
    ///
    /// assert_eq!(inv_a, pse_a); // Nearly equal
    /// ```
    fn pseudo_inv(&self) -> Option<Matrix> {
        let x = self.clone();
        let xt = self.t();

        let xtx = xt % x;
        let inv_temp = xtx.inv();

        match inv_temp {
            None => None,
            Some(m) => Some(m % self.t())
        }
    }
}

// =============================================================================
// Back-end Utils
// =============================================================================

/// Combine separated matrix to one matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix!(1;16;1, 4, 4, Row);
/// let (m1, m2, m3, m4) = a.block();
/// let m = combine(m1,m2,m3,m4);
/// assert_eq!(m, a);
///
/// let b = matrix!(1;16;1, 4, 4, Col);
/// let (n1, n2, n3, n4) = b.block();
/// let n = combine(n1,n2,n3,n4);
/// assert_eq!(n, b);
/// ```
pub fn combine(m1: Matrix, m2: Matrix, m3: Matrix, m4: Matrix) -> Matrix {
    let l = m1.col;
    let c_l = m2.col;
    let r_l = m3.row;

    let r = l + r_l;
    let c = l + c_l;

    let mut m = matrix(vec![0f64; r*c], r, c, m1.shape);

    for idx_row in 0 .. r {
        for idx_col in 0..c {
            match (idx_row, idx_col) {
                (i, j) if (i < l) && (j < l) => {
                    m[(i, j)] = m1[(i, j)];
                },
                (i, j) if (i < l) && (j >= l) => {
                    m[(i, j)] = m2[(i, j - l)];
                },
                (i, j) if (i >= l) && (j < l) => {
                    m[(i, j)] = m3[(i - l, j)];
                },
                (i, j) if (i >= l) && (j >= l) => {
                    m[(i, j)] = m4[(i - l, j - l)];
                },
                _ => (),
            }
        }
    }
    m
}

/// Inverse of Lower matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let a = matrix(c!(1,0,2,1), 2, 2, Row);
/// assert_eq!(inv_l(a), matrix(c!(1,0,-2,1), 2, 2, Row));
///
/// let b = matrix(c!(1,0,0,2,1,0,4,3,1), 3, 3, Row);
/// assert_eq!(inv_l(b), matrix(c!(1,0,0,-2,1,0,2,-3,1), 3, 3, Row));
/// ```
pub fn inv_l(l: Matrix) -> Matrix {
    let mut m = l.clone();

    match l.row {
        1 => l,
        2 => {
            m[(1, 0)] = - m[(1, 0)];
            m
        },
        _ => {
            let (l1, l2, l3, l4) = l.block();

            let m1 = inv_l(l1);
            let m2 = l2;
            let m4 = inv_l(l4);
            let m3 = -((m4.clone() % l3) % m1.clone());

            combine(m1, m2, m3, m4)
        }
    }
}

/// Inverse of upper triangular matrix
///
/// # Examples
/// ```
/// extern crate peroxide;
/// use peroxide::*;
///
/// let u = matrix(c!(2,2,0,1), 2, 2, Row);
/// assert_eq!(inv_u(u), matrix(c!(0.5,-1,0,1), 2, 2, Row));
/// ```
pub fn inv_u(u: Matrix) -> Matrix {
    let mut w = u.clone();

    match u.row {
        1 => {
            w[(0, 0)] = 1f64 / w[(0, 0)];
            w
        },
        2 => {
            let a = w[(0, 0)];
            let b = w[(0, 1)];
            let c = w[(1, 1)];
            let d = a * c;

            w[(0, 0)] = 1f64 / a;
            w[(0, 1)] = - b / d;
            w[(1, 1)] = 1f64 / c;
            w
        },
        _ => {
            let (u1, u2, u3, u4) = u.block();
            let m1 = inv_u(u1);
            let m3 = u3;
            let m4 = inv_u(u4);
            let m2 = -(m1.clone() % u2 % m4.clone());

            combine(m1, m2, m3, m4)
        }
    }
}
