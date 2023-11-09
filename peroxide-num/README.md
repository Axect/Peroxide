# Peroxide-num

`Peroxide-num` is a Rust crate dedicated to providing comprehensive numeric structures and operations, extending generic programming capabilities with a focus on mathematical computations.

## Overview

This crate defines a set of traits for numeric operations, including basic arithmetic, power functions, trigonometric functions, exponential, and logarithmic functions. It is designed to be used with numeric types that implement these traits, allowing for a wide range of mathematical operations.

## Traits

- `PowOps`: Operations related to powers and roots.
- `TrigOps`: Trigonometric functions.
- `ExpLogOps`: Exponential and logarithmic functions.
- `Float`: Define own floating point type (`f32` and `f64` are implemented as default).
- `Numeric`: A comprehensive trait that encompasses all of the above along with standard arithmetic operations.

## Example: Defining a Simple Numeric Type

Below is an example of how you can define your own simple numeric type that implements the `Numeric` trait.

```rust
#[derive(Debug, Clone, Copy, PartialOrd)]
struct SimpleNumber(f64);

impl PowOps for SimpleNumber {
    type Float = Self;

    fn powi(&self, n: i32) -> Self {
        SimpleNumber(self.0.powi(n))
    }

    fn powf(&self, f: Self::Float) -> Self {
        SimpleNumber(self.0.powf(f.0))
    }

    fn pow(&self, f: Self) -> Self {
        SimpleNumber(self.0.powf(f.0))
    }

    fn sqrt(&self) -> Self {
        SimpleNumber(self.0.sqrt())
    }
}

// Implement other required operations for SimpleNumber...
// - Add, Sub, Mul, Div, Neg
// - PowOps (implemented above)
// - TrigOps
// - ExpLogOps

impl Numeric<f64> for SimpleNumber {}
```

This `SimpleNumber` struct wraps a `f64` and implements the `Numeric` trait, making it capable of all the operations defined in the `Peroxide-num` crate.

## Usage

To use this type in your own computations:

```rust
let num = SimpleNumber(2.0);
let result = num.sin(); // Compute the sine of 2.0
println!("{:?}", result); // Should display the sine of 2.0
```

The `Peroxide-num` crate is designed to be flexible and extensible, allowing for easy integration with the larger `Peroxide` ecosystem for scientific computing in Rust.

For more information and advanced usage, please refer to the documentation and examples provided with the crate.
