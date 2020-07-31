//! Stable implementations for nightly-only features
//!
//! # Implemented List
//!
//! * `StableFn` : Make `FnOnce` to stable

/// Stable Fn trait
///
/// # Description
/// Implement `FnOnce` is still nighlty only feature. This trait is alternative to `FnOnce` trait.
pub trait StableFn<T> {
    type Output;
    fn call_stable(&self, target: T) -> Self::Output;
}
