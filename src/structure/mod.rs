//! Main structures for peroxide
//!
//! * Matrix
//! * Vector
//! * Dual
//! * Hyper Dual
//! * Polynomial
//! * DataFrame (only for `dataframe` feature)
//! * Multinomial (not yet implemented)

pub mod dual;
pub mod hyper_dual;
pub mod matrix;
pub mod multinomial;
pub mod polynomial;
pub mod vector;
#[cfg(feature = "dataframe")]
pub mod dataframe;
