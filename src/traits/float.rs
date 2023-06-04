pub trait FloatWithPrecision {
    fn round_with_precision(&self, precision: usize) -> Self;
    fn floor_with_precision(&self, precision: usize) -> Self;
    fn ceil_with_precision(&self, precision: usize) -> Self;
}

impl FloatWithPrecision for f64 {
    fn round_with_precision(&self, precision: usize) -> Self {
        let p = 10f64.powi(precision as i32);
        (self * p).round() / p
    }

    fn floor_with_precision(&self, precision: usize) -> Self {
        let p = 10f64.powi(precision as i32);
        (self * p).floor() / p
    }

    fn ceil_with_precision(&self, precision: usize) -> Self {
        let p = 10f64.powi(precision as i32);
        (self * p).ceil() / p
    }
}

impl FloatWithPrecision for f32 {
    fn round_with_precision(&self, precision: usize) -> Self {
        let p = 10f32.powi(precision as i32);
        (self * p).round() / p
    }

    fn floor_with_precision(&self, precision: usize) -> Self {
        let p = 10f32.powi(precision as i32);
        (self * p).floor() / p
    }

    fn ceil_with_precision(&self, precision: usize) -> Self {
        let p = 10f32.powi(precision as i32);
        (self * p).ceil() / p
    }
}
