use peroxide_num::{PowOps, ExpLogOps, TrigOps, Numeric};
use std::ops::{Neg, Add, Sub, Mul, Div};
use std::f64::consts::PI;

fn main() {
    let a = Vec3D::new(0f64, 0.5 * PI, PI);
    println!("{:?}", a);
    println!("{:?}", -a);
    println!("{:?}", a + a);
    println!("{:?}", a - a);
    println!("{:?}", a * 2f64);
    println!("{:?}", 2f64 * a);
    println!("{:?}", a / 2f64);
    println!("{:?}", 2f64 / a);
    println!("{:?}", a.sin());
    println!("{:?}", a.cos());
    println!("{:?}", a.tan());
    println!("{:?}", a.asin());
    println!("{:?}", a.acos());
    println!("{:?}", a.atan());
    println!("{:?}", a.sinh());
    println!("{:?}", a.cosh());
    println!("{:?}", a.tanh());
    println!("{:?}", a.asinh());
    println!("{:?}", a.acosh());
    println!("{:?}", a.atanh());
    println!("{:?}", a.exp());
    println!("{:?}", a.ln());
}

#[derive(Debug, Copy, Clone)]
struct Vec3D {
    x: f64,
    y: f64,
    z: f64
}

impl Vec3D {
    fn new(x: f64, y: f64, z: f64) -> Vec3D {
        Vec3D { x, y, z }
    }
}

impl Neg for Vec3D {
    type Output = Vec3D;

    fn neg(self) -> Vec3D {
        Vec3D::new(-self.x, -self.y, -self.z)
    }
}

impl Add for Vec3D {
    type Output = Vec3D;

    fn add(self, other: Vec3D) -> Vec3D {
        Vec3D::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl Sub for Vec3D {
    type Output = Vec3D;

    fn sub(self, other: Vec3D) -> Vec3D {
        Vec3D::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Mul for Vec3D {
    type Output = Vec3D;

    fn mul(self, _other: Vec3D) -> Vec3D {
        unimplemented!()
    }
}

impl Div for Vec3D {
    type Output = Vec3D;

    fn div(self, _other: Vec3D) -> Vec3D {
        unimplemented!()
    }
}

impl Add<f64> for Vec3D {
    type Output = Vec3D;

    fn add(self, scalar: f64) -> Vec3D {
        Vec3D::new(self.x + scalar, self.y + scalar, self.z + scalar)
    }
}

impl Sub<f64> for Vec3D {
    type Output = Vec3D;

    fn sub(self, scalar: f64) -> Vec3D {
        Vec3D::new(self.x - scalar, self.y - scalar, self.z - scalar)
    }
}

impl Mul<f64> for Vec3D {
    type Output = Vec3D;

    fn mul(self, scalar: f64) -> Vec3D {
        Vec3D::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Div<f64> for Vec3D {
    type Output = Vec3D;

    fn div(self, scalar: f64) -> Vec3D {
        Vec3D::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }
}

impl Add<Vec3D> for f64 {
    type Output = Vec3D;

    fn add(self, other: Vec3D) -> Vec3D {
        Vec3D::new(self + other.x, self + other.y, self + other.z)
    }
}

impl Sub<Vec3D> for f64 {
    type Output = Vec3D;

    fn sub(self, other: Vec3D) -> Vec3D {
        Vec3D::new(self - other.x, self - other.y, self - other.z)
    }
}

impl Mul<Vec3D> for f64 {
    type Output = Vec3D;

    fn mul(self, other: Vec3D) -> Vec3D {
        Vec3D::new(self * other.x, self * other.y, self * other.z)
    }
}

impl Div<Vec3D> for f64 {
    type Output = Vec3D;

    fn div(self, other: Vec3D) -> Vec3D {
        Vec3D::new(self / other.x, self / other.y, self / other.z)
    }
}

impl PowOps for Vec3D {
    type Float = f64;

    fn pow(&self, _power: Self) -> Self {
        unimplemented!()
    }

    fn powf(&self, power: Self::Float) -> Self {
        Vec3D::new(self.x.powf(power), self.y.powf(power), self.z.powf(power))
    }

    fn powi(&self, power: i32) -> Self {
        Vec3D::new(self.x.powi(power), self.y.powi(power), self.z.powi(power))
    }

    fn sqrt(&self) -> Self {
        Vec3D::new(self.x.sqrt(), self.y.sqrt(), self.z.sqrt())
    }
}

impl TrigOps for Vec3D {
    fn sin(&self) -> Self {
        Vec3D::new(self.x.sin(), self.y.sin(), self.z.sin())
    }

    fn cos(&self) -> Self {
        Vec3D::new(self.x.cos(), self.y.cos(), self.z.cos())
    }

    fn tan(&self) -> Self {
        Vec3D::new(self.x.tan(), self.y.tan(), self.z.tan())
    }

    fn asin(&self) -> Self {
        Vec3D::new(self.x.asin(), self.y.asin(), self.z.asin())
    }

    fn acos(&self) -> Self {
        Vec3D::new(self.x.acos(), self.y.acos(), self.z.acos())
    }

    fn atan(&self) -> Self {
        Vec3D::new(self.x.atan(), self.y.atan(), self.z.atan())
    }

    fn sinh(&self) -> Self {
        Vec3D::new(self.x.sinh(), self.y.sinh(), self.z.sinh())
    }

    fn cosh(&self) -> Self {
        Vec3D::new(self.x.cosh(), self.y.cosh(), self.z.cosh())
    }

    fn tanh(&self) -> Self {
        Vec3D::new(self.x.tanh(), self.y.tanh(), self.z.tanh())
    }

    fn asinh(&self) -> Self {
        Vec3D::new(self.x.asinh(), self.y.asinh(), self.z.asinh())
    }

    fn acosh(&self) -> Self {
        Vec3D::new(self.x.acosh(), self.y.acosh(), self.z.acosh())
    }

    fn atanh(&self) -> Self {
        Vec3D::new(self.x.atanh(), self.y.atanh(), self.z.atanh())
    }

    fn sin_cos(&self) -> (Self, Self) {
        (self.sin(), self.cos())
    }
}

impl ExpLogOps for Vec3D {
    type Float = f64;

    fn exp(&self) -> Self {
        Vec3D::new(self.x.exp(), self.y.exp(), self.z.exp())
    }

    fn ln(&self) -> Self {
        Vec3D::new(self.x.ln(), self.y.ln(), self.z.ln())
    }

    fn log(&self, base: Self::Float) -> Self {
        Vec3D::new(self.x.log(base), self.y.log(base), self.z.log(base))
    }

    fn log2(&self) -> Self {
        Vec3D::new(self.x.log2(), self.y.log2(), self.z.log2())
    }

    fn log10(&self) -> Self {
        Vec3D::new(self.x.log10(), self.y.log10(), self.z.log10())
    }
}

impl Numeric<f64> for Vec3D {}
