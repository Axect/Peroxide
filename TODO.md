# TODO

## 2024.04.14

### Primary

- [ ] Whole new Optimize (Trait based)

- [ ] Pure Rust implementation of Linear Algebra
    - [x] LU (Completely Pivoting)
    - [x] LU (Partial Pivoting)
    - [ ] QR
    - [ ] SVD

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
    - [x] Weighted Uniform
- [ ] Implement special polynomial
    - [x] Legendre
    - [ ] Bessel
    - [ ] Hermite
- [ ] Implement convenient structure of Neural Network
- [ ] Add Statistical regression
    - [ ] Gaussian Kernel
    - [ ] Logistic Kernel
- [ ] Implement more Eigenvalue algorithms
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
- [x] Reduce compile time
  - [x] Replace `proc_macro` for `AD` with ordinary macro or Enum
- [x] Make `csv` optional
- [x] Remove `dual`, `hyperdual` and modify `Real`, `Number` (How to bind `f64` & `AD` effectively?)
- [x] Add more IO options for DataFrame
    - [x] CSV (`csv` feature)
    - [x] NetCDF (`nc` feature)
    - [x] Parquet
- [x] Documentized
    - [x] Vector
    - [x] Matrix
    - [x] Linear Algebra
    - [x] Functional Programming
    - [x] Statistics
    - [x] Interpolation & Spline
    - [x] ODE
    - [x] Macros
    - [x] Optimize
    - [x] Automatic Differentiation
    - [x] DataFrame
- [x] Implement more spline algorithms
- [x] Whole new ODE (trait based)
- [x] Whole new root finding (trait based)
