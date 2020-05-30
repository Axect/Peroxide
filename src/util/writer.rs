//! More convenient matrix writer

pub use self::ToWriter::{Data, Header, Path, Round};
use std::collections::HashMap;
//use std::fs::File;
//use std::io::Write;
//use std::process::exit;
use crate::structure::matrix::Matrix;

#[derive(Debug, Clone, Copy, Hash, PartialOrd, PartialEq, Eq)]
pub enum ToWriter {
    Header,
    Round,
    Data,
    Path,
}

#[derive(Debug, Clone, Copy)]
pub enum Queue {
    Matrix,
    Vector,
}

/// Simple Matrix Writer
///
/// # Necessary Fields
/// * `data: Vec<Matrix>`
/// * `path: String`
///
/// # Option
/// * `header: Vec<String>`
/// * `round: u8`
#[derive(Debug, Clone)]
pub struct SimpleWriter {
    header: Vec<String>,
    round: u8,
    matrices: Vec<Matrix>,
    vectors: Vec<Vec<f64>>,
    path: String,
    queue: Vec<Queue>,
    to_write: HashMap<ToWriter, bool>,
}

impl SimpleWriter {
    pub fn new() -> SimpleWriter {
        let mut default_to_write: HashMap<ToWriter, bool> = HashMap::new();
        default_to_write.insert(Header, false);
        default_to_write.insert(Round, false);
        default_to_write.insert(Data, false);
        default_to_write.insert(Path, false);

        SimpleWriter {
            header: vec![],
            round: 0,
            matrices: vec![],
            vectors: vec![],
            path: "".to_string(),
            queue: vec![],
            to_write: default_to_write,
        }
    }

    pub fn insert_header(&mut self, head: Vec<&str>) -> &mut Self {
        if let Some(x) = self.to_write.get_mut(&Header) {
            *x = true
        }
        self.header = head
            .into_iter()
            .map(|t| t.to_owned())
            .collect::<Vec<String>>();
        self
    }

    pub fn set_round_level(&mut self, nth: u8) -> &mut Self {
        if let Some(x) = self.to_write.get_mut(&Round) {
            *x = true
        }
        self.round = nth;
        self
    }

    pub fn insert_matrix(&mut self, mat: Matrix) -> &mut Self {
        if let Some(x) = self.to_write.get_mut(&Data) {
            *x = true
        }
        self.matrices.push(mat);
        self.queue.push(Queue::Matrix);
        self
    }

    pub fn insert_vector(&mut self, vec: Vec<f64>) -> &mut Self {
        if let Some(x) = self.to_write.get_mut(&Data) {
            *x = true
        }
        self.vectors.push(vec);
        self.queue.push(Queue::Vector);
        self
    }

    pub fn set_path(&mut self, path: &str) -> &mut Self {
        if let Some(x) = self.to_write.get_mut(&Path) {
            *x = true
        }
        self.path = path.to_owned();
        self
    }

    pub fn write_csv(self) {
        unimplemented!()
    }

    // pub fn write_pickle(&self) {
    //     let mut writer: Box<dyn Write>;

    //     // Error handling - Path
    //     if let Some(p) = self.to_write.get(&Path) {
    //         assert!(*p, "No determined path!");
    //     }

    //     // Error handling - Data
    //     if let Some(dat) = self.to_write.get(&Data) {
    //         assert!(*dat, "No inserted data!");
    //     }

    //     match File::create(self.path.clone()) {
    //         Ok(p) => writer = Box::new(p),
    //         Err(e) => {
    //             println!("{:?}", e);
    //             exit(1);
    //         }
    //     }

    //     if let Some(head) = self.to_write.get(&Header) {
    //         if *head {
    //             serde_pickle::to_writer(&mut writer, &self.header, true)
    //                 .expect("Can't write header to pickle");
    //         }
    //     }

    //     let mut queued = self.queue.clone().into_iter();
    //     let mut matrices = self.matrices.clone().into_iter();
    //     let mut vectors = self.vectors.clone().into_iter();

    //     loop {
    //         match queued.next() {
    //             Some(Queue::Matrix) => {
    //                 let mat = matrices.next().unwrap();
    //                 mat.write_pickle(&mut writer)
    //                     .expect("Can't insert matrices");
    //             }
    //             Some(Queue::Vector) => {
    //                 let vec = vectors.next().unwrap();
    //                 vec.write_pickle(&mut writer).expect("Can't insert vectors");
    //             }
    //             None => return,
    //         }
    //     }
    // }
}
