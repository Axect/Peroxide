# TODO

## 2020.12.21

### Primary

- [ ] Pure Rust implementation of Linear Algebra
    - [x] LU (Completely Pivoting)
    - [x] LU (Partial Pivoting)
    - [ ] QR
    - [ ] SVD
- [ ] Add more IO options for DataFrame
    - [x] CSV
    - [x] NetCDF (`nc` feature)
    - [ ] JSON
    - [ ] Arrow IPC
    - [ ] Parquet
- [ ] Remove `dual`, `hyperdual` and modify `Real`, `Number` (How to bind `f64` & `AD` effectively?)
- [ ] Reduce compile time
    - [ ] Replace `proc_macro` for `AD` with ordinary macro or Enum

### Subs

- [ ] Implement `WithJSON` for `DataFrame`
    - [x] `to_json_value`
    - [ ] `from_json_value`
- [ ] Implement various pdf
    - [x] Bernoulli
    - [x] Beta
    - [x] Binomial
    - [ ] Dirichlet
    - [x] Gamma
    - [x] Student's t
    - [x] Uniform
    - [ ] Wishart
- [ ] Implement special polynomial
    - [x] Legendre
    - [ ] Bessel
    - [ ] Hermite
- [ ] Implement convenient structure of Neural Network
- [ ] Documentized
    - [x] Vector
    - [x] Matrix
    - [x] Linear Algebra
    - [x] Functional Programming
    - [x] Statistics
    - [ ] Interpolation & Spline
    - [x] ODE
    - [ ] Macros
    - [ ] Machine Learning
    - [x] Optimize
    - [x] Automatic Differentiation
    - [x] DataFrame
- [ ] Add Statistical regression
    - [ ] Gaussian Kernel
    - [ ] Logistic Kernel
- [ ] Make `csv` optional
- [ ] Make or Use pure Rust plot library
- [ ] Implement more Eigenvalue algorithms
- [ ] Implement more spline algorithms
- [ ] Complex matrix

## Complete

- [x] Can choose API - MATLAB, Python, R
- [x] Implement Plot
- [x] Re-write `numerical` module
- [x] Optimize
    - [x] Linear Regression
    - [x] Non-linear Regression
        - [x] Gauss-Newton (But not yet merged)
        - [x] Gradient Descent
        - [x] Levenberg-Marquardt
- [x] Implement DataFrame
- [x] Implement higher order automatic derivatives
- [x] Generic trait for Automatic differentiation (Create `AD` trait)
- [x] Separate `DataFrame` from `dataframe` feature. (And rename `dataframe` feature to some awesome name)
