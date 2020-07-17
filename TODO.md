# TODO

## 2020.07.17

### Primary

- [ ] Generic trait for Automatic differentiation (Create `AD` trait)
- [ ] Remove `dual`, `hyperdual` and modify `Real`, `Number` (How to bind `f64` & `AD` effectively?)
- [ ] Implement more Eigenvalue algorithms
- [ ] Make or Use pure Rust plot library
- [ ] Separate `DataFrame` from `dataframe` feature. (And rename `dataframe` feature to some awesome name)
- [ ] Make `csv` optional
- [ ] Implement more spline algorithms
- [ ] Complex matrix

### Subs

- [ ] Implement `WithJSON` for `DataFrame`
    - [x] `to_json_value`
    - [ ] `from_json_value`
- [ ] Implement various pdf
    - [x] Bernoulli
    - [x] Beta
    - [ ] Dirichlet
    - [x] Gamma
    - [x] Student's t
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
    - [ ] Automatic Differentiation
- [ ] Add Statistical regression
    - [ ] Gaussian Kernel
    - [ ] Logistic Kernel

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
