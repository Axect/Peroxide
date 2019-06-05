pub trait PowOps {
    type Output;
    fn powi(&self, n: i32) -> Self::Output;
    fn powf(&self, f: f64) -> Self::Output;
    fn sqrt(&self) -> Self::Output;
}

pub trait TrigOps {
    type Output;
    fn sin(&self) -> Self::Output;
    fn cos(&self) -> Self::Output;
    fn tan(&self) -> Self::Output;
    fn asin(&self) -> Self::Output;
    fn acos(&self) -> Self::Output;
    fn atan(&self) -> Self::Output;
    fn sinh(&self) -> Self::Output;
    fn cosh(&self) -> Self::Output;
    fn tanh(&self) -> Self::Output;
    fn asinh(&self) -> Self::Output;
    fn acosh(&self) -> Self::Output;
    fn atanh(&self) -> Self::Output;
    fn sin_cos(&self) -> (Self::Output, Self::Output);
}

pub trait ExpLogOps {
    type Output;
    fn exp(&self) -> Self::Output;
    fn ln(&self) -> Self::Output;
    fn log(&self, base: f64) -> Self::Output;
    fn log2(&self) -> Self::Output;
    fn log10(&self) -> Self::Output;
}

pub trait Real: PowOps + TrigOps + ExpLogOps {
    fn to_f64(&self) -> f64;
    fn from_f64(f: f64) -> Self;
}

// =============================================================================
// Real Traits for f64
// =============================================================================
impl PowOps for f64 {
    type Output = Self;

    fn powi(&self, n: i32) -> Self::Output {
        self.powi(n)
    }

    fn powf(&self, f: f64) -> Self::Output {
        self.powf(f)
    }

    fn sqrt(&self) -> Self::Output {
        self.sqrt()
    }
}

impl TrigOps for f64 {
    type Output = Self;

    fn sin(&self) -> Self::Output {
        self.sin()
    }

    fn cos(&self) -> Self::Output {
        self.cos()
    }

    fn tan(&self) -> Self::Output {
        self.tan()
    }

    fn asin(&self) -> Self::Output {
        self.asin()
    }

    fn acos(&self) -> Self::Output {
        self.acos()
    }

    fn atan(&self) -> Self::Output {
        self.atan()
    }

    fn sinh(&self) -> Self::Output {
        self.sinh()
    }

    fn cosh(&self) -> Self::Output {
        self.cosh()
    }

    fn tanh(&self) -> Self::Output {
        self.tanh()
    }

    fn asinh(&self) -> Self::Output {
        self.asinh()
    }

    fn acosh(&self) -> Self::Output {
        self.acosh()
    }

    fn atanh(&self) -> Self::Output {
        self.atanh()
    }

    fn sin_cos(&self) -> (Self::Output, Self::Output) {
        self.sin_cos()
    }
}

impl ExpLogOps for f64 {
    type Output = Self;

    fn exp(&self) -> Self::Output {
        self.exp()
    }

    fn ln(&self) -> Self::Output {
        self.ln()
    }

    fn log(&self, base: f64) -> Self::Output {
        self.log(base)
    }

    fn log2(&self) -> Self::Output {
        self.log2()
    }

    fn log10(&self) -> Self::Output {
        self.log10()
    }
}

impl Real for f64 {
    fn to_f64(&self) -> f64 {
        *self
    }

    fn from_f64(f: f64) -> Self {
        f
    }
}