//! Statistical Modules
//!
//! * Basic statistical tools - `stat.rs`
//! * Popular distributions - `dist.rs`
//! * Simple Random Number Generator - `rand.rs`
//! * Basic probabilistic operations - `ops.rs`

#[cfg(feature = "special")]
pub mod dist;
pub mod ops;
pub mod rand;
pub mod stat;
