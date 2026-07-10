//! Probability distributions, random sampling, and statistical operations
//!
//! * Basic statistical tools - `stat.rs`
//! * Popular distributions - `dist.rs` (`rand` feature)
//! * Simple Random Number Generator - `rand.rs` (`rand` feature)
//! * Basic probabilistic operations - `ops.rs`

#[cfg(feature = "rand")]
pub mod dist;
pub mod ops;
#[cfg(feature = "rand")]
pub mod rand;
pub mod stat;
