# Peroxide

[![On crates.io](https://img.shields.io/crates/v/peroxide.svg)](https://crates.io/crates/peroxide)
[![On docs.rs](https://docs.rs/peroxide/badge.svg)](https://docs.rs/peroxide/)
[![travis](https://api.travis-ci.org/Axect/Peroxide.svg?branch=master)](https://travis-ci.org/Axect/Peroxide)  
![maintenance](https://img.shields.io/badge/maintenance-actively--developed-brightgreen.svg)

Pure Rust numeric library contains linear algebra, numerical analysis, statistics and machine learning tools with R, MATLAB, Python like macros.

## Latest README version

Corresponding to `0.8.7`.

## Install

* Add next line to your `cargo.toml`
```toml
peroxide = "0.8"
```

## Module Structure

- __src__
  - __bin__ : For test some libraries
    - [dual.rs](src/bin/dual.rs) : Test automatic differentiation
    - [oxide.rs](src/bin/oxide.rs) : Test any
    - [poly.rs](src/bin/poly.rs) : Test polynomial
    - [lotka.rs](src/bin/lotka.rs) : Solve lotka-volterra equation
    - [mlp.rs](src/bin/mlp.rs) : Multi-layer perceptron example
    - [tov.rs](src/bin/tov.rs) : Solve Tolman-Oppenheimer-Volkoff equation
  - [lib.rs](src/lib.rs) : `mod` and `re-export`
  - __ml__ : For machine learning (*Beta*)
    - [mod.rs](src/ml/mod.rs)
    - [reg.rs](src/ml/reg.rs) : Regression tools
  - __macros__ : Macro files
    - [matlab_macro.rs](src/macros/matlab_macro.rs) : MATLAB like macro
    - [mod.rs](src/macros/mod.rs)
    - [r_macro.rs](src/macros/r_macro.rs) : R like macro
  - __numerical__ : To do numerical things
    - [bdf.rs](src/numerical/bdf.rs) : Backward Differentiation Formula
    - [gauss_legendre.rs](src/numerical/gauss_legendre.rs) : Gauss-Legendre 4th order
    - [interp.rs](src/numerical/interp.rs) : Interpolation
    - [mod.rs](src/numerical/mod.rs)
    - [newton.rs](src/numerical/newton.rs) : Newton's Method
    - [ode.rs](src/numerical/ode.rs) : Merge all ODE algorithm to one module
    - [runge_kutta.rs](src/numerical/runge_kutta.rs) : Runge Kutta 4th order
    - [spline.rs](src/numerical/spline.rs) : Natural Spline
    - [utils.rs](src/numerical/utils.rs) : Utils to do numerical things (e.g. jacobian)
  - __operation__ : To define general operations
    - [extra_ops.rs](src/operation/extra_ops.rs)
    - [mod.rs](src/operation/mod.rs)
  - __statistics__ : Statistical Tools
    - [mod.rs](src/statistics/mod.rs)
    - [dist.rs](src/statistics/dist.rs) : Probability distributions
    - [rand.rs](src/statistics/rand.rs) : Wrapper for `rand` crate
    - [stat.rs](src/statistics/stat.rs) : Statistical tools
  - __structure__ : Fundamental data structures
    - [dataframe.rs](src/structure/dataframe.rs) : Not yet implemented
    - [dual.rs](src/structure/dual.rs) : Dual number system for automatic differentiation
    - [matrix.rs](src/structure/matrix.rs) : Matrix
    - [multinomial.rs](src/structure/multinomial.rs) : For multinomial (*Beta*)
    - [mod.rs](src/structure/mod.rs)
    - [polynomial.rs](src/structure/polynomial.rs) : Polynomial
    - [vector.rs](src/structure/vector.rs) : Extra tools for `Vec<f64>`
  - __util__
    - [mod.rs](src/util/mod.rs)
    - [api.rs](src/util/api.rs) : Matrix constructor for various language style 
    - [non_macro.rs](src/util/non_macro.rs) : Primordial version of macros
    - [print.rs](src/util/print.rs) : To print conveniently
    - [useful.rs](src/util/useful.rs) : Useful utils to implement library


## Documentation

There is [Peroxide Gitbook](https://axect.github.io/Peroxide_Gitbook)


## Not yet documentized contents

### Polynomial

```rust
// Peroxide
extern crate peroxide;
use peroxide::*;

fn main() {
    // Declare polynomial
    let a = poly(c!(1,3,2));
    a.print(); // x^2 + 3x + 2
    a.eval(1); // Evaluate when x = 1 -> 6.0
    
    let b = poly(c!(1,2,3,4));       // x^3 + 2x^2 + 3x + 4
    (a.clone() + b.clone()).print(); // x^3 + 3x^2 + 6x + 6
    (a.clone() - b.clone()).print(); // -x^3 - x^2 - 2
    (a.clone() * b.clone()).print(); // x^5 + 5x^4 + 11x^3 + 17x^2 + 18x + 8
    
    let c = poly(c!(1, -1));
    c.print();                       // x - 1
    c.pow(2).print();                // x^2 - 2x + 1
}
```

### Interpolation (Beta)

* Lagrange polynomial interpolation
* Chebyshev nodes

```rust
// Peroxide
extern crate peroxide;
use peroxide::*;

fn main() {
    let a = c!(-1, 0, 1);
    let b = c!(1, 0, 1);

    let l = lagrange_polynomial(a, b);
    l.print(); // x^2
}
```

### Spline (Beta)

* Natural cubic spline

```rust
// Peroxide
extern crate peroxide;
use peroxide::*;

fn main() {
    let x = c!(0.9, 1.3, 1.9, 2.1);
    let y = c!(1.3, 1.5, 1.85, 2.1);

    let s = cubic_spline(x, y);

    for i in 0 .. s.len() {
        s[i].print();
    }
    
    // -0.2347x^3 + 0.6338x^2 - 0.0329x + 0.9873
    // 0.9096x^3 - 3.8292x^2 + 5.7691x - 1.5268
    // -2.2594x^3 + 14.2342x^2 - 28.5513x + 20.2094
}
```

### MATLAB like macro

* `zeros` - zero vector or matrix
* `eye` - identity matrix
* `rand` - random matrix (range from 0 to 1)

### Automatic Differentiation

* Implemented AD with dual number structure.
* Available functions
    * `sin, cos, tan`
    * `pow, powf`
    * `+,-,x,/`
    * `exp, ln`

```rust
extern crate peroxide;
use peroxide::*;

fn main() {
    let a = dual(0, 1); // x at x = 0
    a.sin().print();    // sin(x) at x = 0
    
    // value: 0  // sin(0) = 0
    // slope: 1  // cos(0) = 1
}
```

### Jacobian

* Implemented by AD - Exact Jacobian

```rust
extern crate peroxide;
use peroxide::*;

fn main() {
    let xs = c!(1, 1);
    jacobian(xs, f).print();
    
    //      c[0] c[1]
    // r[0]    6    3
}

// f(t, x) = 3t^2 * x
fn f(xs: Vec<Dual>) -> Vec<Dual> {
    let t = xs[0];
    let x = xs[1];

    vec![t.pow(2) * 3. * x]
}
```

### Ordinary Differential Equation

* Solve 1st order ODE with various methods
* Explicit Method
    * `RK4`: Runge-Kutta 4th order
* Implicit Method
    * `BDF1`: Backward Euler
    * `GL4`: Gauss-Legendre 4th order
    
**Caution**

* input function should have form `(Dual, Vec<Dual>) -> Vec<Dual>`    
    
```rust
// Lotka-Volterra
extern crate peroxide;
use peroxide::*;

fn main() {
    // t = 0, x = 2, y = 1
    let xs = c!(2, 1);
    let rk4_records = solve(lotka_volterra, xs.clone(), (0, 10), 1e-3, RK4);
    let bdf_records = solve(lotka_volterra, xs.clone(), (0, 10), 1e-3, BDF1(1e-15));
    let gl4_records = solve(lotka_volterra, xs, (0, 10), 1e-3, GL4(1e-15));
    //rk4_records.write_with_header("example_data/lotka_rk4.csv", vec!["t", "x", "y"]);
    //bdf_records.write_with_header("example_data/lotka_bdf.csv", vec!["t", "x", "y"]);
    gl4_records.write_with_header("example_data/lotka_gl4.csv", vec!["t", "x", "y"]);
}

fn lotka_volterra(_t: Dual, xs: Vec<Dual>) -> Vec<Dual> {
    let a = 4.;
    let c = 1.;

    let x = xs[0];
    let y = xs[1];

    vec![
        a * (x - x*y),
        -c * (y - x*y)
    ]
}
```

## Version Info

To see [RELEASES.md](./RELEASES.md)

## TODO

To see [TODO.md](./TODO.md)