//! Core data structures: `Matrix`, `Vec<f64>` extensions, `DataFrame`, `Polynomial`, and const-generic forward-mode AD (`Jet<N>`)
//!
//! * Matrix (dense & sparse)
//! * Vector
//! * Automatic differentiation (`Jet<N>`)
//! * Polynomial
//! * DataFrame
//! * Multinomial (not yet implemented)

//pub mod complex;
pub mod ad;
pub mod dataframe;
pub mod matrix;
pub mod multinomial;
pub mod polynomial;
pub mod sparse;
pub mod vector;
