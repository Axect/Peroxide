use std::{ fmt, ops::Index, hash::Hash, collections::HashMap, fmt::Debug };
use structure::matrix::{Matrix, py_matrix};
use std::cmp::{ max, min };
use util::useful::tab;

#[derive(Debug, Clone)]
pub struct DataFrame<T> where T: Hash + Eq + Clone {
    header: HashMap<T, usize>,
    data: Vec<Vec<f64>>,
}

impl<T> fmt::Display for DataFrame<T> where T: Hash + Eq + Clone + Debug {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.spread())
    }
}

#[allow(unused_parens)]
impl<T> DataFrame<T> where T: Hash + Eq + Clone + Debug {
    pub fn new(header: Vec<T>, data: Vec<Vec<f64>>) -> Self {
        let mut hash_header: HashMap<T, usize> = HashMap::new();
        let l = header.len();
        header.into_iter().zip((0 .. l)).for_each(|(t, i)| {
            let result = hash_header.insert(t, i);
            match result {
                None => (),
                Some(_) => panic!("Can't assign value to duplicated key"),
            }
        });
        DataFrame {
            header: hash_header,
            data,
        }
    }

    pub fn find_data(&self, head: T) -> &Vec<f64> {
        let i = self.header.get(&head).unwrap();
        &self.data[*i]
    }

    pub fn to_matrix(&self) -> Matrix {
        py_matrix(self.data.clone())
    }

    pub fn spread(&self) -> String {
        // Find maximum length of data
        let sample = self.data.clone();
        let r: usize = sample.clone()
            .into_iter()
            .map(|x| x.len())
            .fold(0, |x, y| max(x, y));
        let c = self.data.len();

        if r > 100 {
            return format!(
                "Result is too larget! Print this result is not efficient!"
            );
        }        

        let mut space: usize = sample
            .into_iter()
            .map(
                |t| t.into_iter().map(|x| min(format!("{:.4}", x).len(), x.to_string().len())).fold(0, |x,y| max(x,y)), // Choose minimum of approx vs normal
            )
            .fold(0, |x, y| max(x, y))
            + 1;

        if space < 5 {
            space = 5;
        }

        let head_iter = self.header.clone();

        let mut result = String::new();

        result.push_str(&tab("", 5));
        for i in 0 .. c {
            for (head, j) in head_iter.iter() {
                if *j == i {
                    result.push_str(&tab(&format!("{:?}", head), space)); // Header
                    break;
                }
            }
        }
        result.push('\n');

        for i in 0..r {
            result.push_str(&tab(&format!("r[{}]", i), 5));
            for j in 0..c {
                let st1 = format!("{:.4}", self.data[j][i]); // Round at fourth position
                let st2 = self.data[j][i].to_string(); // Normal string
                let mut st = st2.clone();

                // Select more small thing
                if st1.len() < st2.len() {
                    st = st1;
                }

                result.push_str(&tab(&st, space));
            }
            if i == (r - 1) {
                break;
            }
            result.push('\n');
        }

        return result;
    }
}

impl<T> Index<T> for DataFrame<T> where T: Hash + Clone + Eq + Debug {
    type Output = Vec<f64>;

    fn index(&self, index: T) -> &Self::Output {
        self.find_data(index)
    }
}
