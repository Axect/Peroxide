# Peroxide

[![On crates.io](https://img.shields.io/crates/v/peroxide.svg)](https://crates.io/crates/peroxide)
[![On docs.rs](https://docs.rs/peroxide/badge.svg)](https://docs.rs/peroxide/)
[![travis](https://api.travis-ci.org/Axect/Peroxide.svg?branch=master)](https://travis-ci.org/Axect/Peroxide)  
![maintenance](https://img.shields.io/badge/maintenance-actively--developed-brightgreen.svg)

Pure Rust numeric library contains linear algebra, numerical analysis, statistics and machine learning tools with R, MATLAB, Python like macros.

## Latest README version

Corresponding to `0.10.0`.

## Install

* Add next line to your `cargo.toml`
```toml
peroxide = "0.10"
```

## Module Structure

- __src__
  - __bin__ : For test some libraries
  - [lib.rs](src/lib.rs) : `mod` and `re-export`
  - __macros__ : Macro files
    - [matlab_macro.rs](src/macros/matlab_macro.rs) : MATLAB like macro
    - [mod.rs](src/macros/mod.rs)
    - [r_macro.rs](src/macros/r_macro.rs) : R like macro
  - __ml__ : For machine learning (*Beta*)
      - [mod.rs](src/ml/mod.rs)
      - [reg.rs](src/ml/reg.rs) : Regression tools
  - __numerical__ : To do numerical things
    - [bdf.rs](src/numerical/bdf.rs) : Backward Differentiation Formula
    - [gauss_legendre.rs](src/numerical/gauss_legendre.rs) : Gauss-Legendre 4th order
    - [interp.rs](src/numerical/interp.rs) : Interpolation
    - [mod.rs](src/numerical/mod.rs)
    - [newton.rs](src/numerical/newton.rs) : Newton's Method
    - [ode.rs](src/grave/ode.rs) : Merge all ODE algorithm to one module
    - [runge_kutta.rs](src/grave/runge_kutta.rs) : Runge Kutta 4th order
    - [spline.rs](src/numerical/spline.rs) : Natural Spline
    - [utils.rs](src/numerical/utils.rs) : Utils to do numerical things (e.g. jacobian)
  - __operation__ : To define general operations
    - [extra_ops.rs](src/operation/extra_ops.rs)
    - [mod.rs](src/operation/mod.rs)
  - __special__ : Wrapper for `special` crate
    - [mod.rs](src/special/mod.rs)
    - [function.rs](src/special/function.rs) : Special functions
  - __statistics__ : Statistical Tools
    - [mod.rs](src/statistics/mod.rs)
    - [dist.rs](src/statistics/dist.rs) : Probability distributions
    - [ops.rs](src/statistics/ops.rs) : Some probabilistic operations
    - [rand.rs](src/statistics/rand.rs) : Wrapper for `rand` crate
    - [stat.rs](src/statistics/stat.rs) : Statistical tools
  - __structure__ : Fundamental data structures
    - [dataframe.rs](src/structure/dataframe.rs) : Not yet implemented
    - [dual.rs](src/structure/dual.rs) : Dual number system for automatic differentiation
    - [hyper_dual.rs](src/structure/hyper_dual.rs) : Hyper dual number system for automatic differentiation
    - [matrix.rs](src/structure/matrix.rs) : Matrix
    - [multinomial.rs](src/structure/multinomial.rs) : For multinomial (*Beta*)
    - [mod.rs](src/structure/mod.rs)
    - [polynomial.rs](src/structure/polynomial.rs) : Polynomial
    - [vector.rs](src/structure/vector.rs) : Extra tools for `Vec<f64>`
  - __util__
    - [mod.rs](src/util/mod.rs)
    - [api.rs](src/util/api.rs) : Matrix constructor for various language style 
    - [non_macro.rs](src/util/non_macro.rs) : Primordial version of macros
    - [pickle.rs](src/util/pickle.rs) : To handle `pickle` data structure
    - [print.rs](src/util/print.rs) : To print conveniently
    - [useful.rs](src/util/useful.rs) : Useful utils to implement library
    - [writer.rs](src/util/writer.rs) : More convenient write system


## Documentation

There is [Peroxide Gitbook](https://axect.github.io/Peroxide_Gitbook)


## Version Info

To see [RELEASES.md](./RELEASES.md)

## TODO

To see [TODO.md](./TODO.md)