//! Smart pointer of `Vec<f64>`
//!
//! # Usage
//!
//! ```rust
//! extern crate peroxide;
//! use peroxide::*;
//!
//! let a = c!(1, 2, 3, 4);
//! let b = c!(5, 6, 7, 8);
//!
//! let a = a.redox();      // Move Vec to RedoxVector
//! let b = b.redox();      // Move Vec to RedoxVector
//!
//! let c = 2f64 * a - b;   // std::ops are implemented
//! (*c).print();           // Unwrap RedoxVector -> Vec<f64>
//! ```
//!
//! # Implemented Operations
//!
//! * `Add` : With `f64` or `RedoxVector`
//! * `Sub` : With `f64` or `RedoxVector`
//! * `Mul` : With `f64` or `RedoxVector` (`Mul<RedoxVector> for RedoxVector` implies dot product)
//! * `Div` : With `f64`

pub mod redoxable;

