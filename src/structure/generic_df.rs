pub trait DType: Debug + Clone {}

impl DType for usize {}
impl DType for u8 {}
impl DType for u16 {}
impl DType for u32 {}
impl DType for u64 {}
impl DType for isize {}
impl DType for i8 {}
impl DType for i16 {}
impl DType for i32 {}
impl DType for i64 {}
impl DType for f32 {}
impl DType for f64 {}
impl DType for bool {}
impl DType for String {}
impl DType for Char {}

#[derive(Debug, Clone)]
pub trait Column {}

#[derive(Debug, Clone)]
pub struct Series<T: DType> {
    pub name: String,
    pub values: Vec<T>,
}

pub struct DataFrame {
    pub data: Vec<Column>
}

impl<T> Series<T> {
    pub fn new(name: &str, v: Vec<T>) -> Self {
        Self {
            name: name.into(),
            values: v,
        }
    }
}

