# Release 0.39.8 (2025-06-23)

- Implement `LogNormal` distribution
  - `LogNormal(mu: f64, sigma: f64)`
- Fix sampling method for `Gamma`

# Release 0.39.7 (2025-05-27)

- Add some methods for `DataFrame`
  - `filter_by<F: Fn(Scalar) -> bool>(&self, column: &str, f: F) -> anyhow::Result<DataFrame>`
    : Filter rows by a condition on a specific column
  - `mask(&self, mask: &Series) -> anyhow::Result<DataFrame>`
    : Filter rows by a boolean mask
  - `select_rows(&self, indices: &[usize]) -> DataFrame`
    : Select specific rows by indices

# Release 0.39.6 (2025-05-16)

- New ODESolver: `RKF78`
  - Implement `RKF78` method for `ODESolver`

# Release 0.39.5 (2025-04-21)

- New feature `rkyv`
  - Implement `rkyv::{Archive, Serialize, Deserialize}` for `Matrix`, `Polynomial`, `Spline`, `ODE`

# Release 0.39.4 (2025-04-11)

## Optimize `integrate`

- Replace the output signature of `gauss_legendre_table` and `kronrod_table` to `&'static [f64]` to avoid unnecessary allocations.
- Hard code symmetry of weights and nodes into source code to avoid unnecessary allocations.
- New helper function - `compute_gauss_kronrod_sum_stored`
  - Reduce the number of function calls (G+K -> K)
- Change update method of subinterval tolerance (divide by 2 -> divide by sqrt(2))
- These changes improve the performance of `integrate` by 1.2x - 50x (to integrate highly oscillatory functions)

## Update dependencies

- Update `rand` to `0.9`
- Update `rand_distr` to `0.5`

# Release 0.39.3 (2025-03-13)

- Update `puruspe` to `0.4.0`

# Release 0.39.2 (2025-02-06)

- Implement `Broyden` method for `GL4`

# Release 0.39.1 (2025-02-06)

- Add `lambert_w` doc for crate docs [#82](https://github.com/Axect/Peroxide/pull/82) (Thanks to [@JSorngard](https://github.com/JSorngard))

- Add default signature for `linspace!` [#85](https://github.com/Axect/Peroxide/pull/85) (Thanks to [@tarolling](https://github.com/tarolling))

- Fix a bug in `ButcherTableau::step`

- Add another example for ODE (`examples/ode_test_orbit.rs`)

# Release 0.39.0 (2024-11-22) (Yanked)

- Decouple `initial_conditions` from `ODEProblem`
  - Now, we can define `initial_conditions` in solving phase 

# Release 0.38.1 (2024-11-06)

- Fix error in `O3` feature

# Release 0.38.0 (Yanked)

## New features - Complex & Parallel

- `complex` feature
  - Implement complex vector, matrix and integral [#35](https://github.com/Axect/Peroxide/issues/35) (Thanks to [@GComitini](https://github.com/GComitini) and [@soumyasen1809](https://github.com/soumyasen1809))

- `parallel` feature
  - Implement some parallel functions [#72](https://github.com/Axect/Peroxide/issues/72) (Thanks to [@soumyasen1809](https://github.com/soumyasen1809))

## Generic Matrix - MatrixTrait

- Implement `MatrixTrait` for Matrix (Scalar = f64)
- Implement `MatrixTrait` for ComplexMatrix (Scalar = C64)
- `LinearAlgebra` and `solve` depend on `MatrixTrait`

## Other changes

- Update `puruspe` dependency to `0.3.0`, remove `lambert_w` dependency [#79](https://github.com/Axect/Peroxide/pull/79) (Thanks to [@JSorngard](https://github.com/JSorngard))

- Add `hermite_polynomial` and `bessel_polynomial` [#80](https://github.com/Axect/Peroxide/pull/80) (Thanks to [@jgrage](https://github.com/jgrage))

# Release 0.37.9 (2024-07-31)

- Fix inconsistent lambert w function name [#65](https://github.com/Axect/Peroxide/issues/65) (Thanks to [@JSorngard](https://github.com/JSorngard))

# Release 0.37.8 (2024-07-30)

- Integrate with [lambert_w](https://crates.io/crates/lambert_w) crate ([#63](https://github.com/Axect/Peroxide/pull/63)) (Thanks to [@JSorngard](https://github.com/JSorngard))
  - Write flexible wrapper for [lambert_w](https://crates.io/crates/lambert_w)
    ```rust
    pub enum LambertWAccuracyMode {
        Simple,   // Faster, 24 bits of accuracy
        Precise,  // Slower, 50 bits of accuracy
    }

    pub fn lambert_w0(z: f64, mode: LambertWAccuracyMode) -> f64;
    pub fn lambert_wm1(z: f64, mode: LambertWAccuracyMode) -> f64;
    ```

  - Write default Lambert W function for `prelude` (Precise as default)
    ```rust
    use peroxide::prelude::*;

    fn main() {
        lambert_w0(1.0).print(); // Same as fuga::lambert_w0(1.0, LambertWAccuracyMode::Simple)
    }
    ```

# Release 0.37.7 (2024-07-05)

- Bump `pyo3` dependency to `0.22`
- Fix plot functions to be compatible with `pyo3`
- Add B-Spline to README

# Release 0.37.6 (2024-06-19)

## Huge Spline Change

- Generic Spline trait
  - `Spline<T>`: desired output type is `T`
- Split `PolynomialSpline` from `Spline`
  - `CubicSpline` & `CubicHermiteSpline` are now `PolynomialSpline`
  - Implement `Spline<f64>` for `PolynomialSpline`
- Implement B-Spline
  - `BSpline { degree: usize, knots: Vec<f64>, control_points: Vec<Vec<f64>> }`
  - `BSpline::open(degree, knots, control_points)` : Open B-Spline
  - `BSpline::clamped(degree, knots, control_points)` : Clamped B-Spline
- Implement `Spline<(f64, f64)>` for `BSpline`

# Release 0.37.5 (2024-06-10)

- More generic & stable root finding macros (except `Newton`)

# Release 0.37.4 (2024-05-17)

- Public ODE Integrator fields

# Release 0.37.3 (2024-05-01)

- Add Nan/infinite guard to `gauss_kronrod_quadrature` (early exit) ([#59](https://github.com/Axect/Peroxide/pull/59)) (Thanks to [@GComitini](https://github.com/GComitini))
- Add complex feature & complex module ([#35](https://github.com/Axect/Peroxide/issues/35))
- Implement Cubic B-Spline basis functions
  - `UnitCubicBasis`
  - `CubicBSplineBases`

# Release 0.37.2 (2024-04-16)

- Do not include legend box if there is no legend ([#58](https://github.com/Axect/Peroxide/pull/58)) (Thanks to [@GComitini](https://github.com/GComitini))
- Add `rtol` field to `BroydenMethod`
- Implement high-level macros for root finding
  - `bisection!(f, (a,b), max_iter, tol)`
  - `newton!(f, x0, max_iter, tol)` (require `#[ad_function]` attribute)
  - `secant!(f, (a,b), max_iter, tol)`
  - `false_position!(f, (a,b), max_iter, tol)`

# Release 0.37.1 (2024-04-15)

- Implement `BrodenMethod`: Broyden's method (`I>=1, O>=1, T=([f64; I], [f64; I])`)
- Restore citation file

# Release 0.37.0 (2024-04-14)

## Huge Update - Whole new Root finding & anyhow

### Whole new Root finding

- Remove all boilerplates
- Now, `RootFinding` is composed of traits
  - `RootFindingProblem<const I: usize, const O: usize, T>`: Trait for defining and root finding problem
    - `I`: Input dimension
    - `O`: Output dimension
    - `T`: Type of state
  - `RootFinder`: Trait for finding root
    - `BisectionMethod`: Bisection Method (`I=1, O=1, T=(f64, f64)`)
    - `FalsePositionMethod`: False Position Method (`I=1, O=1, T=(f64, f64)`)
    - `NewtonMethod`: Newton Method (`I=1, O=1, T=f64`)
    - `SecantMethod`: Secant Method (`I=1, O=1, T=(f64, f64)`)

### Error handling with anyhow

- Remove `thiserror` dependency
- Add `anyhow` for error handling
- Change error handling in `ODE`, `Spline`, `WeightedUniform`

# Release 0.36.4 (2024-04-11)

- More generic Butcher tableau
  - Now, you can use `ButcherTableau` for non-embedded methods too
- More ODE integrators
  - `RALS3, RALS4, RK5, BS23`

# Release 0.36.3 (2024-04-10)

- Hotfix : Fix `GL4` algorithm

# Release 0.36.2 (2024-04-10)

- Now, you can report current states if your constraints are violated.
  - `ODEError::ConstraintViolation` -> `ODEError::ConstraintViolation(f64, Vec<f64>, Vec<f64>)`
  - for detailed information, see [docs for ODEError](https://axect.github.io/Peroxide_Doc/peroxide/numerical/ode/enum.ODEError.html)
- Add docs for `ODEError`

# Release 0.36.1 (2024-04-09)

- Fix all warnings in peroxide
- Change redundant method
  - `Vec<f64>::resize` -> `Vec<f64>::reshape`
- Error handling for concatenation
  - `cbind` & `rbind` now returns `Result<Matrix, ConcatenateError>`
- New non-macro utils
  - `column_stack(&[Vec<f64>]) -> Result<Matrix, ConcatenateError>`
  - `row_stack(&[Vec<f64>]) -> Result<Matrix, ConcatenateError>`
  - `rand_with_rng(usize, usize, &mut Rng) -> Matrix`
- Generic Butcher tableau trait (now for embedded Runge-Kutta methods)

  ```rust
  pub trait ButcherTableau {
      const C: &'static [f64];
      const A: &'static [&'static [f64]];
      const BH: &'static [f64];
      const BL: &'static [f64];
  
      fn tol(&self) -> f64;
      fn safety_factor(&self) -> f64;
      fn max_step_size(&self) -> f64;
      fn min_step_size(&self) -> f64;
      fn max_step_iter(&self) -> usize;
  }
  ```

  - Implement `ODEIntegrator` for `ButcherTableau`
    - Just declare `ButcherTableau` then `step` is free

  - Three available embedded Runge-Kutta methods
    - `RKF45`: Runge-Kutta-Fehlberg 4/5th order
    - `DP45`: Dormand-Prince 4/5th order
    - `TSIT45`: Tsitouras 4/5th order

# Release 0.36.0 (2024-04-08)

## Huge Update - Error handling & Whole new ODE

### Error handling

- Add `thiserror` for error handling
- Implement errors for cubic spline & cubic hermite spline.
- Implement errors for weighted uniform distribution & PRS.

### Seedable sampling

- Now, all distribution has `sample_with_rng` method.
- There are two wrappers for `SeedableRng`
  - `smallrng_from_seed` : Performant but not secure
  - `stdrng_from_seed` : Performant enough and secure enough

### Whole new ODE

- Remove all boilerplates.
- Now, `ODE` is composed of traits.
  - `ODEProblem`: Trait for defining and ODE problem.
  - `ODEIntegrator`: Trait for integrating ODE.
    - `RK4`: Runge-Kutta 4th order
    - `RKF45`: Runge-Kutta-Fehlberg 4/5th order
    - `GL4`: Gauss-Legendre 4th order
    - You can implement your own integrator.
  - `ODESolver`: Trait for solving ODE.
    - `BasicODESolver`: Basic ODE solver - define range of t, initial step size and integrate it.
    - You can implement your own solver.
- For more information, see [docs for ode](https://axect.github.io/Peroxide_Doc/peroxide/numerical/ode/index.html).
    
# Release 0.35.1 (2024-03-29)

- Add `PlotType` for `Plot2D`
  - `PlotType::Scatter`
  - `PlotType::Line` (default)
  - `PlotType::Bar`

# Release 0.35.0 (2024-03-29)

## Change some plot functions

- Now you can set marker, line style, color, alpha option for specific element.
  - `set_marker(vec![(usize, Marker)])` : `usize` is index of element (image or pair)
  - `set_line_style(vec![(usize, LineStyle)])`
  - `set_color(vec![(usize, String)])`
  - `set_alpha(vec![(usize, f64)])`

# Release 0.34.7 (2024-03-11)

## More updates for `plot` feature

- Make legend optional (Now, no legend is available)
- Implement `set_line_style`. Here are available line styles.
  - `LineStyle::Solid`
  - `LineStyle::Dashed`
  - `LineStyle::Dotted`
  - `LineStyle::DashDot`
- Implement `set_color`
- Implement `set_alpha`
- More markers.

## Getter for ODE

- Add explicit getter for `ExplicitODE` and `ImplicitODE` for various fields.

# Release 0.34.6 (2024-03-01)

## Big updates for `plot` feature

* Add `auto-initialize` flag for `pyo3`
* Add `scienceplots` support. Here are available styles.
  * `PlotStyle::Default` : default matplotlib style - no scienceplots required
  * `PlotStyle::Science` : scienceplots default style - scienceplots required
  * `PlotStyle::Nature` : nature style - scienceplots required
  * `PlotStyle::IEEE` : IEEE style - scienceplots required
* Implement `xscale, yscale, xlim, ylim` for `Plot2D`
* You can check these features in [Peroxide Gallery](https://github.com/Axect/Peroxide_Gallery/tree/master/Plot/plot_feature)

# Release 0.34.5 (2024-02-08)

* Derive `Serialize` and `Deserialize` for `CubicHermiteSpline`

# Release 0.34.4 (2024-01-28)

* Derive `Serialize` and `Deserialize` for `Matrix`
* Remove explicit implementation for `Default` for `Shape`

# Release 0.34.3 (2023-11-25)

* Update `peroxide-num` to `v0.1.4`
* Implement `ExpLogOps, PowOps, TrigOps` and `Numeric<f64>` for `Matrix`

# Release 0.34.2 (2023-11-20)

## New sub-crate : `peroxide-num`

* Add new sub-crate : `peroxide-num`
* Change all dependencies of `ExpLogOps, PowOps, TrigOps` to `peroxide-num`

## Fix errata

* R example in `structure/matrix` (#56) (Thanks to [@rdavis120](https://github.com/rdavis120))

## Fix old syntax

* Fix old syntax - e.g. explicit `into_iter`, `Vec::with_capacity` & `set_len`

# Release 0.34.1 (2023-08-03)

## Modify Statistics of `WeightedUniform`

* Modify `self.sum` to compatible with definition [Weighted Uniform Distribution](https://axect.github.io/posts/006_prs/#22-weighted-uniform-distribution)
* Modify `mean` & `var` to compatible with definition

## Add more `Vec<f64>` algorithms

* Implement `arg_min`, `max`, `min`

# Release 0.34.0 (2023-06-13)

## Change Gauss-Kronrod quadrature 

* Adapt max iteration number to Gauss-Kronrod quadrature
    * Arguments of all methods related with Gauss-Kronrod quadrature are changed.
    * e.g. `G7K15(1e-15)` -> `G7K15(1e-15, 20)` (20 is maximum iteration)

* Example
    ```rust
    use peroxide::fuga::*;

    fn main() {
        let f_integral = integrate(f, (0f64, 1f64), G7K15R(1e-4, 20));
        f_integral.print();
    }

    fn f(x: f64) -> f64 {
        x.powi(2)
    }
    ```

# Release 0.33.5 (2023-06-13)

* Implement Gauss-Kronrod quarature with relative error
    * `G7K15R`, `G10K21R`, `G15K31R`, `G20K41R`, `G25K51R`, `G30K61R` can be used.
* Reduce warning messages

# Release 0.33.4 (2023-06-04)

* Implement `Statistics` for `WeightedUniform` ([#55](https://github.com/Axect/Peroxide/pull/55)) (Thanks to [@samnaughtonb](https://github.com/samnaughtonb))
* New trait: `FloatWithPrecision`
    * `fn round_with_precision(&self, precision: usize) -> Self`
    * `fn floor_with_precision(&self, precision: usize) -> Self`
    * `fn ceil_with_precision(&self, precision: usize) -> Self`
* New utils:
    * `fn seq_with_precision(start, end, step, precision: usize) -> Vec<f64>`
    * `fn linspace_with_precision(start, end, length, precision: usize) -> Vec<f64>`

# Release 0.33.2 (2023-04-03)

* Implement necessary traits for `ConfusionMatrix`
    * `#[derive(Debug, Clone, PartialEq)]`
* Bump up dependencies version
    * `netcdf` : `0.7.0` -> `0.8.1`
    * `arrow` : `0.14` -> `0.17.0`
    * `pyo3`: `0.17` -> `0.18`

# Release 0.33.1 (2023-03-14)

* Implement `ConfusionMatrix` in `statistics::stat`
    * Implement all metrics in [wikipedia](https://en.wikipedia.org/wiki/Confusion_matrix)

# Release 0.33.0 (2023-03-08)

* Delete `build.rs` to remove any explicit linkages to specific BLAS implementations ([#54](https://github.com/Axect/Peroxide/issues/54)) (Thanks to [@gfaster](https://github.com/gfaster))

# Release 0.32.1 (2022-11-04)

* Make an option for choosing compression method for parquet
    * At `fuga` : `fn write_parquet(&self, path: &str, compression: CompressionOptions)`
    * At `prelude` : `fn write_parquet(&self, path:&str)` (Default: `CompressionOptions::Uncompressed`)

# Release 0.32.0 (2022-11-03)

## DataFrame meets Parquet

* Add `parquet` feature
* Add `WithParquet` trait and implement it for `DataFrame`
    * `fn write_parquet(&self, path: &str) -> Result<(), Box<dyn Error>>`
    * `fn read_parquet(path: &str) -> Result<Self, Box<dyn Error>>`
    * Update `DataFrame` docs

# Release 0.31.8 (2022-10-11)

* Change debug procedure for stop condition of `ODE` ([#52](https://github.com/Axect/Peroxide/issues/52)) (Thanks to [@tchamelot](https://github.com/tchamelot))
    * Add `fn has_stopped(&self) -> bool` for `ODE` struct

# Release 0.31.7 (2022-09-21)

* Fix bug in `linspace` ([#51](https://github.com/Axect/Peroxide/issues/51))
* Change print scheme of `Vec<float>`
    * Now, floating number in `Vec` is printed by `fmt_lower_exp(2)`

# Release 0.31.6 (2022-07-15)

* Add `*_with_cond` for `Spline` trait (See [Truncated Cubic - Peroxide Gallery](https://github.com/Axect/Peroxide_Gallery/tree/master/Numeric/truncated_cubic) for an example)
    * `eval_with_cond<F: Fn(f64) -> f64>(&self, x: f64, cond: F) -> f64` : Evaluate with custom condition
    * `eval_vec_with_cond<F: Fn(f64) -> f64 + Copy>(&self, v: [&f64], cond: F) -> f64` : Evaluate vector with custom condition

# Release 0.31.5 (2022-07-06)

* New trait - `LowerExpWithPlus`, `UpperExpWithPlus`
    * Now, we can print `132.45` as `1.3245e+2` via `132.45.fmt_lower_exp(4)`
    * Now, we can print `132.45` as `1.3245E+2` via `132.45.fmt_upper_exp(4)`
* Change print scheme of `DataFrame`
    * Now, floating number in DataFrame is printed by `fmt_lower_exp(2)`

# Release 0.31.4 (2022-06-15)

* Fix bug in `rref` ([#50](https://github.com/Axect/Peroxide/pull/50))

# Release 0.31.3 (2022-05-29)

* Fix bug in `linspace` and `logspace`
    * Now `linspace(0, 0, 1)` returns `[0f64]` instead of `[NaN]`
    * Now `logspace(0, 0, 1, 10)` returns `[1f64]` instead of `[NaN]`

# Release 0.31.2 (2022-05-29) (Yanked)

* Fix assertion of `util::non_macro::seq`
* Implement numpy like `logspace`

# Release 0.31.1 (2022-05-24)

* Fix a bug in `spread` of DataFrame

# Release 0.31.0 (2022-05-19)

## Significant changes

### Splines

* New trait `Spline`
    * Move `CubicSpline::eval` to `Spline::eval`
    * Move `CubicSpline::polynomial` to `Spline::polynomial_at`
    * Move `CubicSpline::number_of_polynomials` to `Spline::number_of_polynomials`
    * Add `Spline::eval_vec`
    * Add `Spline::get_ranged_polynomials`
* Implement Cubic Hermite spline
    * Add struct `CubicHermiteSpline`
    * Implement slope estimation algorithms
        * `SlopeMethod::Akima`
        * `SlopeMethod::Quadratic`
* Modify `CubicSpline` **(Important!)** 
    * Change argument type
        * `from_nodes(node_x: &[f64], node_y: &[f64]) -> Self`
    * (For developer) Remove `CubicSpline::ranged` (use `util::useful::zip_range` instead)
* Add docs for `numeric/spline.rs`

### Calculus

* Rename `Calculus::diff` to `Calculus::derivative` **(Important!)**

## Minor updates

* Add `util::useful::{gen_range, zip_range}`
* Add `structure::poly::Calculus::integrate`

# Release 0.30.15 (2022-05-02)

* Update `puruspe` to `0.2.0` (Fix a bug in gamma function)

# Release 0.30.14 (2022-03-22)

* New distribution : `WeightedUniform`
* New sampling method : `PRS` (Piecewise Rejection Sampling)
* Documentation for these new methods will be added later.

# Release 0.30.13 (2022-02-21)

* Control `lambda` of `LevenbergMarquardt` (#49)
    * Add `set_lambda_init` & `set_lambda_max` (#49)

# Release 0.30.12 (2022-02-21)

* More flexible root finding
    * Add getter methods for `RootFinder`
    * Add getter methods for `RootState`

# Release 0.30.11 (2022-02-10)

* Implement Cholesky decomposition in `O3` feature.
* Implement symmetricity check method - `is_symmetric` for Matrix. 

# Release 0.30.10 (2022-02-05)

* Update `netcdf` dependency to `0.7`
    * Fix `nc` feature issue - not compatible with hdf5 version 1.12.0
* Update `pyo3` dependency to `1.15`
* Update `float-cmp` dev dependency to `0.9`

# Release 0.30.9 (2021-05-26)

* Add more trigonometric ops
    * `asin_acos(&self) -> (Self, Self)`
    * `asinh_acosh(&self) -> (Self, Self)`
* Update dependencies
    * `blas` : `0.21.0` -> `0.22.0`
    * `lapack` : `0.17.0` -> `0.19.0`

# Release 0.30.8 (2021-05-21)

* Fix errata in `col_map`, `row_map`

# Release 0.30.7 (2021-05-14)

* Change signature of `cubic_spline`
    * Originally, `(Vec<f64>, Vec<f64>) -> Vec<Polynomial>`
    * Now, `(&Vec<f64>, &Vec<f64>) -> CubicSpline` 
* Add Truncated SVD
    * Add `truncated(&self)` method for `SVD`

# Release 0.30.6 (2021-03-31)

* Fix a bug in quantile of `statistics/stat.rs`

# Release 0.30.5 (2021-03-26)

* Update docs
    * `prelude/mod.rs` : Update default numerical integration method
* Update `matrixmultiply` dependency
    * Add `threading` feature
    * Enhance matrix multiplication performance : See [matmul](https://github.com/Axect/Scientific_Bench/tree/master/Linear_Algebra/matmul)

# Release 0.30.4 (2021-03-01)

* Update dependencies
    * `rand` : 0.7 -> 0.8
    * `rand_distr` : 0.3 -> 0.4
    * `matrixmultiply` : 0.2 -> 0.3
    * `netcdf` : 0.5 -> 0.6
    * `blas` : 0.20 -> 0.21
    * `lapack` : 0.16 -> 0.17
    * `pyo3` : 0.12 -> 0.13

# Release 0.30.3 (2021-01-26)

* Automatic generated gradient & hessian via `proc_macro`
    * Currently only support `Fn(f64) -> f64`

```rust
use peroxide::fuga::*;

fn main() {
    f(2f64).print();        // x^3     = 8
    f_grad(2f64).print();   // 3 * x^2 = 12
    f_hess(2f64).print();   // 6 * x   = 12
}

#[ad_function]              // generates f_grad, f_hess
fn f(x: f64) -> f64 {
    x.powi(3)               // x^3
}
```

# Release 0.30.2 (2021-01-21)

* Implement Gauss-Kronrod Quadrature
    * G7K15
    * G10K21
    * G15K31
    * G20K41
    * G25K51
    * G30K61
* Now, prelude's default integration is `G7K15(1e-16)`

# Release 0.30.1 (2021-01-20)

* Implement Chebyshev polynomial
* Implement more higher order Gauss-Legendre Quadrature (Up to 30)

# Release 0.30.0 (2021-01-14)

## Independence day from Dual

* Replace all `Dual`, `HyperDual`, `Number` with `AD`
    * Remove `Dual, HyperDual, Number`
* It affected all of numerical functions
    * `numerical/root.rs`
    * `numerical/ode.rs`
    * `numerical/optimize.rs`
    * `numerical/utils.rs`
* Also many traits are changed
    * `traits/num.rs`
    * `traits/pointer.rs`

## Minor Changes

* No more default in `util::non_macro::concat`
* Now, `VecOps` has default implmentations - It required `FPVector`

# Release 0.29.1 (2021-01-13)

* Implements all numerical operations of `AD`
    * Inverse trigonometric: `asin, acos, atan`
    * Inverse hyperbolic: `asinh, acosh, atanh`
    * Power of AD: `pow(&self, other: AD) -> Self`

# Release 0.29.0 (2021-01-13)

## Diet

* Remove all `proc_macro` (Remove dependency of `peroxide-ad`)
  * Now, `AD` are enums
    * `AD0(f64)`
    * `AD1(f64, f64)`
    * `AD2(f64, f64, f64)`
* `csv` becomes optional
  * Remove dependency of `csv, serde, iota, ...`

# Release 0.28.2 (2021-01-12)

* Reduce compile time via `watt` integration
    * Now, `peroxide-ad` is pre-compiled to wasm

# Release 0.28.1 (2020-12-22)

* Fix dimension error of `apply` in `O3` feature
* Import `structure/dataframe.rs` into `prelude`
* Update version of `README.md`

# Release 0.28.0 (2020-12-21)

## Rebirth of DataFrame

* Now, `DataFrame` can contain multiple type columns. Refer to [dataframe](https://peroxide.surge.sh/structure/dataframe/index.html).
* `DataFrame` is merged default feature. No more `dataframe` feature required. Thus, `dataframe` feature is removed.
* But if you want to `netcdf` file format, then `nc` feature is required.

## Minor changes

* Fix errata in `README.md`
* Remove unnecessary imports

# Release 0.27.1 (2020-10-28)

* Add doc for `Matrix::qr` & `Matrix::svd`
* Enhance doc for `Matrix::pseudo_inv`
* Implement `pseudo_inv` via `svd` (`O3` feature only)

# Release 0.27.0 (2020-10-12)

* Update dependencies
    * `indexmap` : 1.5 -> 1.6
    * `pyo3` : 0.11 -> 0.12
* Add `print` for `&Vec<T>`

# Release 0.26.4 (2020-09-24)

* Fix dimension error of `SVD.u`

# Release 0.26.3 (2020-09-24)

## SVD

* Add `svd` to `Matrix` (Only available in `O3` feature)

## Minor changes

* Remove `packed_simd` dependency (Fix `packed_simd` error)
* Add `O3` version of `qr` (using `lapack_dgeqrf`)

# Release 0.26.2 (2020-09-23)

* Add more sugar for `Vec<f64>`
    * `ConvToMat` : Convert `Vec<f64>`
    * `to_col` : To Column matrix
    * `to_row` : To Row matrix

# Release 0.26.1 (2020-09-23)

* Add new methods for `Matrix`
    * `submat(&self, start: (usize, usize), end: (usize, usize))` : Return submatrix
    * `subs_mat(&mut self, start, end, &Matrix)` : Substitute submatrix

# Release 0.26.0 (2020-08-28)

* Update `netcdf` dependencies
    * Now, use `netcdf = 0.5`
* Add new methods for `DataFrame`
    * `head_print(&self, n: usize)` : Return n lines from head
    * `tail_print(&self, n: usize)` : Return n lines before tail

# Release 0.25.8 (2020-08-28)

* Change licenses : BSD-3-Clause -> MIT OR Apache-2.0
* Update version of dependencies

# Release 0.25.7 (2020-08-08)

* Implement more Vector Products
    * `cross`
    * `outer`
* Implement more Matrix Products
    * `kronecker`
    * `hadamard`

# Release 0.25.6 (2020-08-06)

* Add assertion for matrix multiplications

# Release 0.25.5 (2020-08-04)

* Implement Binomial distribution

# Release 0.25.4 (2020-08-01)

* Reduce compile time
    * Reduce order of `AD{i}` : 10 -> 5

# Release 0.25.3 (2020-07-31)

## Root Finding is available

* Add `numerical/root.rs` (See [docs](https://peroxide.surge.sh))
    * Low-level API
    * High-level API
* Add `ADLift<F, T>` for lifting generic `AD` function

# Release 0.25.2 (2020-07-23)

* Impl `std::ops` with `f64` for `AD{i}` and vice versa
* Increase Accuracy of Spline Extension (Thanks to [schrieveslaach](https://github.com/schrieveslaach))

# Release 0.25.1 (2020-07-17)

* Integrate `src/structure/ad.rs` into `prelude`
* Add generic automatic differenitation trait - `AD`

# Release 0.25.0 (2020-07-15)

## Higher order Automatic Differentiation

* Add `structure/ad.rs`
* Add `peroxide-ad` (`proc_macro` for AD)
* Implement `AD1` ~ `AD10` (Upto 10th order)
* Modify `traits/num.rs`
    * Change some methods to provided methods

# Release 0.24.5 (2020-06-18)

* Add `get_env(&self)` in `ExplicitODE`
* Add `get_env(&self)` in `ImplicitODE`

# Release 0.24.4 (2020-06-18)

* Fix `set_header` error of `DataFrame`

# Release 0.24.3 (2020-06-18)

## Solve ODE with Environment

* Add `Environment` trait in `numerica/ode.rs`
    * `ODE` -> `ODE<E: Environment>`
    * `ExplicitODE` -> `ExplicitODE<E: Environment>`
        * `f: Fn(&mut State<f64>)` -> `f: Fn(&mut State<f64>, &E)`
        * Add `set_env(E)`
    * `ImplicitODE` -> `ImplicitODE<E: Environment>`
        * `f: Fn(&mut State<Dual>)` -> `f: Fn(&mut State<Dual>, &E)`
        * Add `set_env(E)`

# Release 0.24.2 (2020-06-18) (Yanked)

# Release 0.24.1 (2020-06-16)

* Fetch `prelude` with new Linear algebra
    * Add `SimpleLinearAlgebra`

# Release 0.24.0 (2020-06-16)

## Whole new LU

* Fix error in `LinearAlgebra::lu`
    * Add `gepp, gecp` for partial pivoting and complete pivoting
    * Peroxide chooses `gecp` default
* No more `unwrap`
    * `lu` returns `PQLU` directly
    * `inv` returns `Matrix` directly
    * `pseudo_inv` returns `Matrix` directly

## Now, we can solve!

* Implement two solve algorithms
    * LU decomposition via Gaussian elimination with Complete pivoting (Stable)
    * WAZ decomposition (Unstable)

### Example

```rust
#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = ml_matrix("1 2;3 4");
    let b = c!(3, 7);
    
    a.solve(&b, LU).print();    // [1, 1]
    a.solve(&b, WAZ).print();   // [1, 1]
}
```

# Release 0.23.2 (2020-06-10)

* Fix errata in `Sub` of `Redox`

# Release 0.23.1 (2020-06-10)

## More Tests

* Add `tests/linalg.rs` : It compares `peroxide` and `julia` with `test_data/*.nc`

## More sugar

* Add `traits/sugar.rs`
    * `VecOps` : Vector operation with vectors and scalars
    * `Scalable` : Easy to resize vector or matrix and also concatenation

## More functional

* Add `col_reduce`, `row_reduce`

# Release 0.23.0 (2020-06-01)

**[Caution!] Huge Update!**

## Change module re-exporting

No more direct re-exporting. Below code is not allowed.
```rust
extern crate peroxide;
use peroxide::*;
```  

Now, peroxide has two re-export options - `prelude` and `fuga`.

* `prelude`: To use simple
* `fuga`: To control numerical algorithms

For example,
```rust
// Prelude
#[macro_use]
extern crate peroxide;
use peroxide::prelude::*;

fn main() {
    let a = c!(1, 2, 3);
    assert_eq!(a.norm(), 14f64.sqrt());
}
```

```rust
// Fuga
#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = c!(1, 2, 3);
    assert_eq!(a.norm(Norm::L2), 14f64.sqrt());
}
```

## New module: `traits`

* Remove `operation`
* Create `traits`

`traits` contains below submodules.

* `fp.rs` : Functional Programming toolbox
* `general.rs` : General algorithms
* `math.rs` : Mathematical traits
* `mutable.rs` : Mutable toolbox
* `num.rs` : `Real` & `Number` & Various `Ops`
* `pointer.rs` : `Redox<T: Vector>` & `MatrixPtr`

## Specify Rust edition

* From Ver `0.23.0`, peroxide uses 2018 edition.

## Change `Matrix` - `Vec<f64>` multiplication

* From ver 0.23.0, `Matrix * Vec<f64> = Vec<f64>` and vice versa.

## Minor changes

* Replace `norm_l*()` as `norm(Norm::L*)`
* Move `interp::lagrange_polynomial`, `special::legendre::legendre_polynomial` to `structure/polynomial.rs`
* Remove `LinearOps` in `structure/matrix.rs`
    * Replace `to_matrix(&self)` with `Into<Matrix>`
    * Replace `from_matrix(&self)` with `Into<Vec<f64>>`
    * Move `transpose(), t()` to `Matrix::transpose(), Matrix::t()`
* Remove `VecOps` in `structure/vec.rs`
    * Replace `add, sub, s_mul` with `traits::math::Vector`
    * Replace `norm()` with `traits::math::Normed`
    * Replace `dot()` with `traits::math::InnerProduct`
    * Use `Redox<T: Vector>` rather than `VecOps` (Refer to `traits/pointer.rs`)
* Remove `special/legendre.rs`
* Add `gemv`, `gevm` in `structure/matrix.rs`

# Release 0.22.0 (2020-05-25)

## Add Numerical integration

* `numerical/integral.rs`
    * Newton Cotes quadrature - `integrate(f, (a, b), NewtonCotes(usize))`
    * Gauss Legendre quadrature - `integrate(f, (a, b), GaussLegendre(usize))`

## More `MutOps`

* `swap_with_perm` : Swap with permutation

## More `CubicSpline` (By [schrieveslaach](https://github.com/schrieveslaach))

* `polynomial(&self, x: T) -> Polynomial` : Returns a reference the `Polynomial` at the given point `x`.

## More `Vector::norm` (By [nateckert](https://github.com/Nateckert))

* Add more `norm` for `Vec<f64>`: `norm_l1(&self), norm_l2(&self), norm_linf(&self), norm_lp(&self)`

## TODO

* Gaussian elimination with LU decomposition
* `Perm * Matrix` : Syntactic sugar for `swap_with_perm(_, Row)`
* `Matrix * Perm` : Syntactic sugar for `swap_with_perm(_, Col)`
* More numerical integrations
* Unify vector norms
* Make `csv` to `optional`

# Release 0.21.7 (2020-04-22)

## More `LinearAlgebra`

* QR Decomposition
    * Add `qr(&self) -> QR` in `LinearAlgebra` trait
    * Add `QR` to represent QR decomposition
* Reduced Row Echelon Form
    * Add `rref(&self) -> Matrix` in `LinearAlgebra` trait

## More Efficient

* Modify Polynomial evaluate algorithm via Horner's Method (Thanks to [Nateckert](https://github.com/Nateckert))

# Release 0.21.6 (2020-04-11)

* Create `util/wrapper.rs` : Wrapper for other crates.
    * Trait `SampleRNG` : Extract random sample from `Vec<T>`
        * `sample(&self, n: usize) -> Vec<Self::Item>`
* More `Printable`
    * `Vec<usize>`, `Vec<u32>`, `Vec<u64>`
    * `Vec<isize>`, `Vec<i32>`, `Vec<i64>`

# Release 0.21.5 (2020-04-02)

* Fix a bug for `det`, `inv` & `pseudo_inv` for large matrix
    * Set precision of `lu` to more lower : `1e-7` to `1e-40`
    
## TODO

* QR decomposition
* Effective pseudo inverse algorithm using QR
* Enhance performance of `Matrix`

# Release 0.21.4 (2020-03-21)

## Important

* New dependency - `matrixmultiply`
* Change default matrix multiplication behavior - depend on `matrixmultiply`
* If size of matrix is smaller than `1000 x 1000`, `default` is faster than `O3`
* New function - `gemm` : Wrapper of `dgemm` of `matrixmultiply`
    * `gemm(alpha, A, B, beta, C)` : `C = alpha * A * B + beta * C`  

# Release 0.21.3 (2020-03-21)

* Add `operation/row_ops.rs`
    * Add `RawMatrix`
    * `row_ptr(&self, usize) -> Vec<*const f64>`
    * `col_ptr(&self, usize) -> Vec<*const f64>`
* Add `as_slice, as_mut_slice` for `Matrix`
* Add `ptr_to_vec` in `util/low_level`

# Release 0.21.2 (2020-03-17)

## Important

* Add `Eigen`
    * Implement *jacobi method*

## Example
```rust
extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = MATLAB::new("1 2; 2 3");
    let eigen = eigen(&a, Jacobi);
    let (eig_val, eig_vec) = eigen.extract();

    eig_val.print();
    eig_vec.print();
}
```

# Release 0.21.1 (2020-03-03)

## Important

* Modify `PowOps`
    * Rename `powf(&self, Self)` to `pow(&self, Self)`
    * Create `powf(&self, f64)`
    
## Other enhancement

* Implement `std::ops` for `RedoxVector`
* Implement `PowOps` for `RedoxVector`
* Update Documents

# Release 0.21.0 (2020-03-01)

## Independence day of `special`

## Important

* Remove dependencies - `special`, `special-fun`
    * Now, use special functions & distributions for WASM.
* New dependency - [puruspe](https://github.com/Axect/puruspe) (**PUR**e **RUS**t **SPE**cial function library)
* Re-implement special functions by `puruspe`
    - [x] `ln_gamma`
    - [x] `gamma`
    - [x] `inc_gamma`
    - [x] `inv_inc_gamma`
    - [x] `beta`
    - [x] `inc_beta`
    - [x] `inv_inc_beta`
    - [x] `erf`
    - [x] `erfc`
    - [x] `inv_erf`
    - [x] `inv_erfc`

# Release 0.20.3 (2020-02-28)

* Add `Div` for `Matrix`

# Release 0.20.2 (2020-02-27)

## Hotfix
* Revert `cfg(feature="special")` as `cfg(feature="specials")` for WASM.
    * `specials = ["special", "special-fun"]`
    
## TODO
* Implement special functions by pure Rust.

# Release 0.20.1 (2020-02-26)

* Add Accessoires to `CubicSpline` (By [schrieveslaach](https://github.com/schrieveslaach))
    - return number of polynomials
    - access element through `std::ops::Index`
* Move `docs` domain : [https://peroxide.surge.sh](https://peroxide.surge.sh)
* Add more R macros
    * `rt`
    * `dt`
    * `pt`

# Release 0.20.0 (2020-02-25)

## Important

* New dependencies
    * `order-stat` : For ordered statistics
    * `float-cmp` : To compare floating numbers conveniently
    * `special-fun` : To use additional special functions
* `OrderedStat`
    * Now, we can calculate quantile (perfectly matched with R quantile)
    * Implemented list
        - [x] Type1
        - [x] Type2
        - [ ] Type3, ... , Type 9
* Remove `special` feature - Now, `special`, `special-fun` are necessary dependencies.
* New method over `RNG` trait. - `cdf`

## Other Enhancement

### Added

* Some additional R macros
    * `rnorm`
    * `dnorm`
    * `prorm`
* Some additional special functions
    * `inc_beta` : Regularized incomplete Beta integral
    * `inc_gamma` : Regularized incomplete lower Gamma integral
    * `hyp2f1` : Hypergeometric function

### Modified

* Update `REAME.md`
* Fix comments of `plot.rs`

# Release 0.19.4 (2020-01-30)

* Modify `vector`
    * Replace last `native` with `O3`
    * Add `sum` to `VecOps`
* Modify documentation of `dist`
    * Remove checkbox
    * Replace `~` with `<del>` tag

# Release 0.19.3 (2019-12-21)

* Remove travis-ci error in `examples`
* Remove some legacy codes
* No differences with `0.19.2`

# Release 0.19.2 (2019-12-21)

* Remove test errors of `dataframe.rs`

# Release 0.19.1 (2019-12-21) (Yanked)

* New dependency - `rand_distr`
    * Now, any distributions are depend on `rand_distr`
* Add `StudentT` distribution
* Rename `example` to `examples`
    * Now use `cargo build --examples`. But it should require all features so, not recommended.
    Instead of this, you can see real examples in [Peroxide Gallery](https://github.com/Axect/Peroxide_Gallary).


# Release 0.19.0 (2019-12-19)

* New dependency in `dataframe` feature - `json`
* Add `WithJSON` trait in `dataframe.rs`
* Implement `WithJSON` for `DataFrame`
    - [x] `to_json_value`
    - [ ] `from_json_value`

# Release 0.18.6 (2019-12-05)

* Add missing trait bound for `Real` trait
* Fix error of `plot.rs`
    * Change python library dependency : `pylab` -> `matplotlib.pyplot`

# Release 0.18.5 (2019-11-30)

* **[Hotfix]** Fix spacing bug for `r>100` cases.

# Release 0.18.4 (2019-11-29)

* Fix spacing constraint - space should be larger than length of key

# Release 0.18.3 (2019-11-29)

* Fix global spacing of `dataframe` to column-wise spacing.

# Release 0.18.2 (2019-10-23)

* Improve generic of `jacobian`

# Release 0.18.1 (2019-10-07)

* Fix limitation of key length: 5 -> unlimited

# Release 0.18.0 (2019-10-06)

## Important

* Rename feature : `oxidize` -> `O3`

## Other enhancement

### Added

* Cubic spline structure (By [schrieveslaach](https://github.com/schrieveslaach))
* Add documentation for `DataFrame`

### Modified

* `structure::DataFrame::WithNetCDF`
    * Modify `read_nc(file_path: &str, header: Vec<&str>)` -> `read_nc(file_path: &str)`
    * Add `read_nc_by_header(file_path: &str, header: Vec<&str>)`

# Release 0.17.3 (2019-10-03)

* Implement `Statistics` for `DataFrame`

# Release 0.17.2 (2019-10-02)

* With surge, new documentation web sites available - [peroxide.info](http://peroxide.info)

# Release 0.17.0 (2019-10-01)

## DataFrame

### Added

* `structure/dataframe.rs`

### New dependency

* `indexmap`
* `netcdf`

### New feature

* `dataframe`

### Implement level for DataFrame

- [x] Pretty print
- [x] Able to print different length
- [x] Convert to Col matrix
- [x] Initialize with header
- [x] Call by header (`Index` with header)
- [ ] Call by row
- [x] Insert pair
- [x] Insert row
- [x] `IndexMut` with header
- [ ] Any column operations
- [x] Read from csv
- [x] Write to csv
- [x] Read from netcdf
- [x] Write to netcdf

# Release 0.16.5 (2019-09-26)

* Remove `NumError` for `Number`
* Add `PartialOrd` for `Number, Dual, HyperDual`

# Release 0.16.4 (2019-09-26)

## Safe Optimization

### Modified 

#### Remove unnecessary assertions

* `structure::dual::ExpLogOps::ln`
* `structure::hyperdual::Div`
* `structure::hyperdual::ExpLogOps::ln`

#### Safe optimization

* `numerical::optimize`
    * `F: Fn(&Vec<f64>, Vec<Number>) -> Option<Vec<Number>>`

# Release 0.16.3 (2019-09-23)

## More mutable operations

### Added

* `util::low_level`
    * `unsafe fn copy_vec_ptr(dst: &mut Vec<*mut f64>, src: &Vec<f64>)`
    * `unsafe fn swap_vec_ptr(lhs: &mut Vec<*mut f64>, rhs: &mut Vec<*mut f64>)`

### Modified

* `structure::matrix::Matrix::swap` -> `operation::mut_ops::MutMatrix::swap`
    * Now `swap` becomes mutable function

# Release 0.16.2 (2019-09-21)

* Now, we can insert pair of data to `plot`
    * `fn insert_pair(&self, pait: (Vec<f64>, Vec<f64>)) -> &mut Self`

# Release 0.16.1 (2019-09-19)

* More generic `optimize`
    * function pointer -> `Box<F>`

# Release 0.16.0 (2019-09-14)

## Now, peroxide meets low level

### Added
* `operation::mut_ops::MutMatrix`
    * `fn col_mut(&mut self, idx: usize) -> Vec<*mut f64>`
    * `fn row_mut(&mut self, idx: usize) -> Vec<*mut f64>`
* `structure::Matrix::FP::{col_mut_map, row_mut_map}`
    
### Modified
* `structure::matrix::Matrix`
    * `fn ptr(&self) -> *const f64`
    * `fn mut_ptr(&self) -> *mut f64`
    * `Index for Matrix`
    * `IndexMut for Matrix`

### Example
```rust
fn main() {
    // ===================================
    // Low Level
    // ===================================
    let mut a = ml_matrix("1 2; 3 4");
    a.print();
    //      c[0] c[1]
    // r[0]    1    2
    // r[1]    3    4

    unsafe {
        let mut p: Vec<*mut f64> = a.col_mut(1); // Mutable second column
        for i in 0 .. p.len() {
            *p[i] = i as f64;
        }
    }

    a.print();
    //      c[0] c[1]
    // r[0]    1    0
    // r[1]    3    1
    
    // ===================================
    // High Level
    // ===================================
    let mut b = ml_matrix("1 2 3; 4 5 6");
    b.col_mut_map(|x| x.normalize());
    b.print();
    //        c[0]   c[1]   c[2]
    // r[0] 0.2425 0.3714 0.4472
    // r[1] 0.9701 0.9285 0.8944
}
```

# Release 0.15.3 (2019-09-12)

* Smart pointer of Vector - `RedoxVector`
* SIMD integrated - `packed_simd`
* Revive Gauss-Legendre 4th order ODE solver

# Release 0.15.2 (2019-09-04)

* Come back to Travis-CI (only default feature)
* Legendre Polynomial (`legendre_polynomial(n: usize) -> Polynomial`)
* Gauss-Legendre Quadrature based on numerical table (up to `n=16`)

# Release 0.15.1 (2019-08-29)

* More optimization of `vector.rs`
* Can use `openblas` feature in `ode.rs`

# Release 0.15.0 (2019-08-27)

* **[Important]** More features
    * You can choose next features
        * `openblas` : BLAS & LAPACK back-end
        * `plot` : Plot with matplotlib (depends on `pyo3` - should use nightly compiler)
* More BLAS
    * Vector & Vector Addition/Subtraction (`daxpy`)
    * Vector & Scalar multiplication (`dscal`)
    * Vector & Vector dot product (`ddot`)
    * Vector Euclidean norm (`dnrm2`)
* More LAPACK
    * QR decomposition (`dgeqrf`)
        * get Q,R (`dorgqr`)
    * get Condition number (`dgecon`)

# Release 0.14.0 (2019-08-24)

* **[Important]** Finally, BLAS integrated.
    * You can choose `blas, lapack` by `native` features.
        ```sh
        cargo build --features native
        ```
    * Default features are no BLAS (Pure Rust)
        ```sh
        cargo build
        ```
    * To use BLAS in Rust can be bothered. Should refer to [OpenBLAS for Rust](https://github.com/Axect/Issues/tree/master/Rust)

    * BLAS implemented Ops
        - [x] Matrix & Matrix Addition/Subtraction (`daxpy`)
        - [x] Matrix & Scalar Addition/Subtraction (`daxpy`)
        - [x] Matrix & Matrix Multiplication (`dgemm`)
        - [x] Matrix & Vector Multiplication (`dgemv`)
    * LAPACK implemented ops
        - [x] LU Factorization (`dgetrf`)
        - [x] Inverse by LU (`dgetri`)
        - [x] Solve by LU (`dgetrs`)
* Move unnecessary `bin` files to example directory

# Release 0.13.0 (2019-08-11)

* **[Important]** Remove inefficient things
    * Remove `util/pickle.rs`
    * Remove `serde`, `serde-pickle` dependencies
* Optimize original code
    * Matrix
        * `change_shape`
        * `col(usize), row(usize)`
        * `diag`
        * `subs_col, subs_row`
        * `Add<Matrix>`
        * `Mul<Matrix>` - Block matrix multiplication

# Release 0.12.4 (2019-08-09)

* Add `normalize` for `Vec<f64>`
* Add `col_map` for `Matrix`
* Add `row_map` for `Matrix`

# Release 0.12.3 (2019-07-31)

* Add `get_error` for `Optimizer`
    * Now, we can see root mean square error of `Optimizer`
* Optimizer documents
* Change non-reasonable syntax
    * `set_legends` -> `set_legend` (in `util/plot.rs`)

# Release 0.12.2 (2019-07-30)

* Add `set_markers` for `Plot2D`
    * Can support `Point`, `Circle`, `Line`

# Release 0.12.1 (2019-07-30)

* Fix error of `powf` in `PowOps` for `Number`
    * Cover all cases
* Change output type of `optimize`
    * Original: `Matrix`
    * Changed: `Vec<f64>`

# Release 0.12.0 (2019-07-29)

* **[Important]** Change the definition of `powf` in `PowOps`
    * `powf(&self, Self) -> Self`
    * Apply this definition to
        * `Dual`
        * `HyperDual`
        * `Number`
    * And remove `PowOps` of `Vec<f64>`

# Release 0.11.8 (2019-07-29) 

* More utils
    * `max<T>(v: Vec<T>) -> T where T: PartialOrd + Copy + Clone`
    * `min<T>(v: Vec<T>) -> T where T: PartialOrd + Copy + Clone`
* Implement Levenberg-Marquardt Algorithm (Refer to `bin/optimize.rs`)
* `to_diag` for `Matrix`
    * `to_diag` : Extract diagonal matrix from a matrix
* Add non-linear regression algorithms (Optimizations) to `numerical/optimize.rs`
    * Gradient Descent
    * Gauss Newton (Not yet implemented)
    * Levenberg Marquardt
* Add julia-like macros
    * `hstack!` : Vectors to Column matrix
    * `vstack!` : Vectors to Row matrix

# Release 0.11.7 (2019-07-24)

* Modify `jacobian`
    * Receive `Fn(Vec<Number>) -> Vec<Number>`
* Remove `bdf.rs`, `gauss_legendre.rs`
* More extended functional programming for `Vec`
    * Extend `zip_with` to any real vector

# Release 0.11.6 (2019-07-22)

* Add `extract` for `PQLU`
    * `extract: PQLU -> (p, q, l, u)`

# Release 0.11.5 (2019-07-20)

* Update whole documentations - part 2.

# Release 0.11.4 (2019-07-17)

* For convenience
    * New method - `to_vec` of `Matrix`
        * `Matrix` to `Vec<Vec<f64>>`
    * Change input type of `set_legends` of `plot.rs`
        * `Vec<String>` -> `Vec<&str>`
* Update whole documentations - part 1.

# Release 0.11.3 (2019-07-14)

* Fix performance issue in `0.11.2`
    * Fix non-efficient parts of `take` of `Matrix`

# Release 0.11.2 (2019-07-14)

* Exclude `bin` directory
* Apply `stop_condition` for `ODE`

# Release 0.11.1 (2019-07-13)

* Remove dependency of `inline-python`

# Release 0.11.0 (2019-07-13)

* **[Important]** Now, only nightly support (Because of `pyo3`)
* Integrate with `inline-python`
* Update `README.md`

# Release 0.10.6 (2019-07-12)

* Now, we can draw & save plot
    * New dependency - `pyo3`
    * Plot is implemented in `util/plot.rs`

# Release 0.10.5 (2019-07-09)

* Modify `integrate` of `numerical/ode.rs`
    * Now, `integrate` provides [`param` | `values`] matrix.
* Change `rand` crate version - `0.7.0`

# Release 0.10.4 (2019-07-01)

* Add `Number` enum for generic function input to ODE.
    * `Number` is composed of `F(f64), D(Dual), E(NumError)`
    * It replace `Real` trait

# Release 0.10.3 (2019-06-27)

* Modify `State<T>` structure. (`state` -> `value`)
    ```rust
    pub struct State<T> {
        param: T,
        value: Vec<T>,
        deriv: Vec<T>
    }
    ```
* Modify `ExplicitODE`
    * Remove `count` field
* Gitbook & README update

# Release 0.10.2 (2019-06-27)

* **[Important!]** Re-define ode structure
    * Great improve UI - Like `SimpleWriter`
    * Implement `ExMethod::Euler`
    * Implement `ExMethod::RK4`
* Fix bug in `SimpleWriter`
    * Fix "always header" bug

# Release 0.10.1 (2019-06-25)

* `std::ops` for `&Dual`
* `std::ops` for `&HyperDual`

# Release 0.10.0 (2019-06-24)

## Huge Updates for Ver 0.10.0

* **[Important!]** Remove `Rem` for `Matrix` (Thanks to [russellb23](https://github.com/russellb23))
* **[Important!]** Change `Mul` for `Matrix`
    * `Mul<Matrix> for Matrix` is Matrix multiplication!
* Now, we can use `std::ops` for `&Matrix`
    ```rust
    extern crate peroxide;
    use peroxide::fuga::*;

    fn main() {
        let a = ml_matrix("1 2;3 4");

        (&a + &a).print();
        (&a * &a).print();
        (&a - &a).print();
    }
    ```
* `Real` Trait is appeared & implemented for some types.
    * `Real` for `f64`
    * `Real` for `Dual`
    * `Real` for `HyperDual`
    ```rust
    extern crate peroxide;
    use peroxide::*;
    
    fn main() {
        let x_f64 = 2f64;
        let x_dual = dual(2, 1);
        let x_hyper = hyper_dual(2, 1, 0);
    
        f(x_f64).print();
        f(x_dual).print();
        f(x_hyper).print();
    }
    
    fn f<T: Real>(x: T) -> T {
        return x.powi(2);
    }
    ```
* Implement dot product for `Vec<Dual>` (Thanks to [russellb23](https://github.com/russellb23))

# Release 0.9.4 (2019-05-13)

* Add `log(&self, base: f64)` to `ExpLogOps`

# Release 0.9.3 (2019-05-13)

* Implement `Printable` to `OPDist<T>, TPDist<T>`
* New trait `ParametricDist` for `OPDist<T>, TPDist<T>` (Just extract parameters)
* New module: `util/writer.rs`
    * You can write pickle file with pipelines.

# Release 0.9.2 (2019-05-03)

* Implement Arnoldi iteration & Gram-schmidt (Not yet merged)
    * `bin/arnoldi.rs`
    * `bin/schmidt.rs`
* Add `Debug, Clone` to `OPDist<T>, TPDist<T>`

# Release 0.9.1 (2019-04-15)

* Add `zeros_shape`, `eye_shape` to `util/non_macro.rs`
* Fix `Matrix::from_index`
    * You should use index function which returns `f64`

# Release 0.9.0 (2019-04-08)

* Modify `Pickle` trait - Allow multiple data to one pickle file
    * `write_single_pickle` : Just write vector or matrix to one pickle file.
    * `write_pickle(&self, writer: &mut Write)`
        ```rust
        extern crate peroxide;
        use peroxide::*;
        use std::fs::File;
        use std::io::Write;

        fn main () {
            let mut w: Box<Write>; 
            match File::create(path) {
                Ok(p) => writer = Box::new(p),
                Err(e) => (),
            }

            let a = ml_matrix("1 2;3 4");
            a.write_pickle(&mut w).expect("Can't write pickle file");
        }
        ```
* Implement matrix norm (usage: `a.norm(<input_norm>)`)
    * `PQ(p, q)` : `L_pq` norm
    * `One` : `L_1` norm
    * `Infinity`: `L_∞` norm
    * `Frobenius` : Frobenius norm (= `PQ(2,2)`)

# Release 0.8.13 (2019-04-01)

* `HyperDual` for 2nd order Automatic Differentiation

# Release 0.8.12 (2019-03-31)

* Implement Tri-Diagonal Matrix Algorithm
    * Add `tdma` to `numerical/utils.rs`

# Release 0.8.11 (2019-03-31)

* Modify `matrix/lu` to correct doolittle algorithm

# Release 0.8.10 (2019-03-25)

* Add two dependencies in `Cargo.toml`
    * `serde`
    * `serde_pickle`
* Add `pickle.rs` to `util`
    * Write `Vec<f64>` to pickle file easily
    * Write `Matrix` to pickle file (Caution: It depends on shape of matrix)
* Fix all warnings from compiler
    
# Release 0.8.9 (2019-03-24)

* New constructor in `matrix.rs`
    * Matrix from index operations - `from_index<F>(F, (usize, usize))`
* Update print of Matrix
    * Extend limit from `10x10` to `100x10` or `10x100` or `20x20`
* Fix bug of `take` of `FPMatrix`
    * Early return if size is smaller than row or column

# Release 0.8.8 (2019-03-23)

* Add `ops.rs` to `statistics`
    * Factorial: `factorial(n)`
    * Permutation: `P(n,r)`
    * Combination: `C(n,r)`
    * Combination with Repetition: `H(n,r)`

# Release 0.8.7 (2019-03-08)

* Add constraint to uniform distribution
* Reduce & modify `README.md`
    * Add missing modules to module structure
    * Remove `Usage` section (Move to [Peroxide Gitbook](https://axect.github.io/Peroxide_Gitbook))
    

# Release 0.8.6 (2019-03-05)

* Add [Peroxide Gitbook](https://axect.github.io/Peroxide_Gitbook) link to `README`
* Fix `statistics/rand.rs`, `statistics/dist.rs`
    * pub use rand crate to private use
    * Now you can use `rand(usize, usize)` function in `util/non_macro.rs`

# Release 0.8.5 (2019-02-17)

* Fix bugs of `cbind`, `rbind`
* Add Linear Discriminant (Least Square) example

# Release 0.8.4 (2019-02-14)

* Fix complete pivoting
* Add `det` test in `tests`
* Add `tov` in `bin` - Tolman-Oppenheimer-Volkoff equation

# Release 0.8.3 (2019-02-13)

* Fix error of `Sub<Dual> for f64`
* Fix error of `Div<Dual> for f64`

# Release 0.8.2 (2019-02-12)

* Bump `rand` dependency to 0.6 (Thanks to [koute](https://github.com/koute))
* Fix error of `powf` operation of `dual`
    * Now, it works fine. 

# Release 0.8.1 (2019-02-04)

* Fix errors of test

# Release 0.8.0 (2019-02-04) (Yanked)

* Fix `write`, `write_with_header`
    * Move `round` parameter to `write_round`, `write_with_header_round`
* Add `solve_with_condition` to `ode.rs`
    * Now, you can give stop condition to ode solver.

# Release 0.7.7 (2019-01-27)

* Add various distributions in `dist.rs`
    * `Bernoulli(mu)`
    * `Beta(a, b)`

# Release 0.7.6 (2019-01-21)

* Modify `write`, `write_with_header`
    * Now there is round option

# Release 0.7.5 (2019-01-21)

* Fix error of `bdf.rs`

# Release 0.7.4 (2019-01-21)

* Modify `bdf.rs`
    * Put `max_iter = 10`
    * Simplify non-autonomous jacobian
    
# Release 0.7.3 (2019-01-19)

* Move distributions(`Uniform`, `Normal`) from `rand.rs` to `dist.rs`
    * Now `Uniform` & `Normal` are enums
    * Remove `Uniform::new` & `Normal::new`
* Add `special/function.rs`
    * Add `gaussian`
    

# Release 0.7.2 (2019-01-18)

* Implement `GL4` - Gauss-Legendre 4th order
    * Add `GL4(f64)` to `ODEMethod`

# Release 0.7.1 (2019-01-17)

* Add `take` to `FP` trait for Matrix
* Add `skip` to `FP` trait for Matrix
* Fix `fmt::Display` of `Matrix`
    * If larger than 10x10 -> Only print 10x10 part
* Add `ode.rs` to `numeric`
    * Add `solve` - numerical solve ODE
    * Now you can choose two methods
        * `RK4`
        * `BDF1`
* Change `rk4`
    * All functions should have form
        * `f(Dual, Vec<Dual>) -> Vec<Dual>`
* Fix error of `spread`

# Release 0.7.0 (2019-01-15)

* Modify matrix declaration
    * `p_matrix` -> `py_matrix`
    * `m_matrix` -> `ml_matrix`
    * Add `r_matrix` (same as `matrix`)
* Add `util/api.rs`
    * Can choose various coding style
        * `MATLAB`
        * `PYTHON`
        * `R`
* Remove `CreateMatrix`
    * Deprecated `Matrix::new` -> Use `matrix` instead

# Release 0.6.15 (2019-01-12)

* Update `matrix.rs`
    * Add `p_matrix`, `m_matrix`
        * Pythonic matrix
        * MATLAB matrix
    * Add `write_with_header` for `matrix.rs`
        * Now, can write matrix with header

# Release 0.6.14 (2019-01-05)

* Add `runge_kutta.rs`
    * Implement RK4 algorithm for Non-autonomous equation

# Release 0.6.13 (2019-01-03)

* Add `grave`
* Move `rok4a.rs` to `grave`

# Release 0.6.12 (2019-01-03)

* Fix error of `Div` for `Dual`

# Release 0.6.11 (2018-12-31)

* Add `rok4a.rs`
    * Now, deprecated
* Add `non_auto_jacobian` to utils
    * TODO: Generalize jacobian
* Add `bdf.rs`
    * Implement Backward Euler Method

# Release 0.6.10 (2018-12-27)

* Add comfortable tools for `Vec<Dual>`
* Add `jacobian` in `numerical/utils`
* Add `newton` in `numerical`
    * Newton-Raphson Method 

# Release 0.6.9 (2018-12-26)

* **Fix error of `inv`**
    * Reverse order of permutations

# Release 0.6.8 (2018-12-25)

* Update `Dual`
    * Also add `Ops<Dual> for f64`

# Release 0.6.7 (2018-12-24)

* Add `multinomial.rs`
    * Implement `print`, `eval`
    * TODO: Partial eval?
* Update `Dual`
    * Add `Add<f64> for Dual`
    * Add `Sub<f64> for Dual`
    * Add `Mul<f64> for Dual`
    * Add `Div<f64> for Dual`

# Release 0.6.6 (2018-11-30)

* Update `FPVector`
    * Add `filter, take, drop`

# Release 0.6.5 (2018-11-29)

* Update `read`
    * Move `read` to `Matrix::read`
    * Can set `delimiter`

# Release 0.6.4 (2018-11-28)

* Add `pseudo_inv` method for `Matrix`
* New `useful.rs` in util
    * Move `tab, quot_rem, nearly_eq` from `matrix.rs` to `useful.rs`
    * Move `choose_*` from `polynomial.rs` to `useful.rs` 
* Fix error of `VectorOps` - `dot`

# Release 0.6.3 (2018-11-28)

* Fix typo of `fmt::Display` for `Polynomial` 
* Fix module structures - Thanks to md-file-tree

# Release 0.6.2 (2018-11-26)

* Implement Horner's Algorithm - `Div` for Polynomial
* TODO: Fix `print` of Polynomial (Fixed!)
* Add `Dual` for Automatic Differentiation to structure
    * structure
        * matrix
        * vector
        * polynomial
        * dual
* Make `operation` directory & add `extra_ops.rs`

# Release 0.6.1 (2018-11-18)

* Fix `README`

# Release 0.6.0 (2018-11-18)

* Implement `Calculus` for `Polynomial`
* Re-construct all module structures
    * structure
        * matrix.rs
        * vector.rs
        * polynomial.rs
    * statistics
        * stat.rs
        * rand.rs
    * macros
        * r_macro.rs
        * matlab_macro.rs
    * util
        * print.rs
* Add `numerical` directory
    * Add interp.rs
        * Lagrange Polynomial (`lagrange_polynomial`)
        * Chebyshev Nodes (`chebyshev_nodes`)
    * Add spline.rs
        * Natural Cubic Spline (`cubic_spline`)
* Impl `pow` for Polynomial
* Fixed `fmt::Display` for Polynomial

# Release 0.5.8 (2018-11-16)

* Add `print.rs` for print any values conveniently
    * Implement `print` for Vector
    * Implement `print` for Matrix
    * Implement `print` for `f32, f64, usize, u32, u64, i32, i64`
* Add `poly.rs` to deal `Polynomial`
    * Implement `fmt::Display` for Polynomial
    * Add `new`
    * Implement `eval` for Polynomial
    * Implement `Neg, Add, Sub, Mul<T>` for Polynomial
    * Implement `Mul, Add<T>, Sub<T>, Div<T>` for Polynomial

# Release 0.5.7 (2018-11-13)

* Change gaussian generate method
    * `marsaglia_polar` to `ziggurat`
* Add comments and examples to `rand.rs`

# Release 0.5.6 (2018-11-13)

* Add `linspace` to `matlab_macro`
* Fixed `linspace` export error
* Add `rand.rs`
    * Generic `Rand` structure
    * `sample` method
    * Marsaglia Polar
    * `Rand` to `Uniform` and `Normal`

# Release 0.5.5 (2018-11-06)

* Extend `matrix` macro to single valued matrix
* Make `lm`
* And also make `lm` macro - `lm!(y ~ x)`
* Make `LinearOps` Trait - But not necessary

# Release 0.5.4 (2018-11-05)

* Add badges to README
* Fix README - add cargo.toml
* Modify `std::ops` for Matrix
    * `f64` to generic
    * Add comments
    * Matmul for `Matrix` vs `Vector` vice versa

# Release 0.5.3 (2018-11-05)

* Add `eye` to `matlab_macro`
* Extend `zeros` to matrix
* Fix `cov` for `Vec<f64>` - not consume anymore
* Add `cor`
* Update `README`

# Release 0.5.2 (2018-11-03)

* Add `matlab_macro`

# Release 0.5.1 (2018-11-02)

* Add `read`
    * Can read matrix from csv
* Add comment to `write`, `read`
* Fix all README

# Release 0.5.0 (2018-11-01)

* Add `write` for `Matrix`
    * Can write matrix to csv!

# Release 0.4.9 (2018-10-30)

* Modify `block`, `inv_u`, `combine`
    * Just change code syntax

# Release 0.4.8 (2018-10-30)

* Modify `lu`
    * Just change code syntax

# Release 0.4.7 (2018-10-30)

* Add `IndexMut` implementation for `Matrix`
* Modify `Rem`
    * Just using `IndexMut`
    * Very fast!

# Release 0.4.6 (2018-10-29)

* Fixed `block` & `combine`
    * Only squared matrices -> Every matrices

# Release 0.4.5 (2018-10-26)

* More add R-like macro
    * `cbind`
    * `rbind`
* README update

# Release 0.4.4 (2018-10-26)

* Refactor structure
    * Move all macro to `r_macro.rs`

# Release 0.4.3 (2018-10-24)

* Add `stat.rs`
    * `mean, var, sd`
* Modify `spread`
    * Fix bugs of all cases
    * Use modern syntax

# Release 0.4.2 (2018-10-24)

* Fix `Cargo.toml`

# Release 0.4.1 (2018-10-24)

* Replace `lu` with `plu`
* Make Algorithm trait for Vector
    * `rank, sign, arg_max`
* Change `PartialEq` for `Matrix`
    * Add `nearly_eq` function
    * Use `nearly_eq` for `Matrix::eq`
* Add `swap`
* Make `PQLU` structure
* Remove `pivot`, `to_perm`
* Replace `plu` with `lu`
    * `lu` returns `Option<PQLU>`
* Enhance error handling with `lu, det, inv`
* Complete Pivoting LU Decomposition
* Fix error of `lu` - initialize `u`

# Release 0.4.0 (2018-10-23)

* Remove non-necessary comments
* Remove `vec2mat, mat2vec`
* Change `col`, `row` functions
    * `col, row` returns `Vec<f64>`
* Add `diag`
* Add `det`
* Add `reduce` in `vector_macro`
* Add `inv_l`, `inv_u`
* Add `block`, `combine`
* Fix error of `block`, `combine`
* Fix error of `inv_l`
* Add `inv`

# Release 0.3.1 (2018-10-21)

* Remove `Vector` struct
    * Replace with `vector_macro`
    * `c!` & `seq!`
* Make R-like matrix macro
    * `matrix!(1;4;1, 2, 2, Row)`

# Release 0.3.0 (2018-10-20)

* Vector
* `seq` : moved from matrix to vector
* Rename `Generic` trait - `CreateMatrix`

# Release 0.2.5 (2018-10-19)

* LU Decomposition

# Release 0.2.4 (2018-10-19)

* `matrix` function - Same as R
* Fix `README.md`
* More documentation

# Release 0.2.3 (2018-10-18)

* `seq` function - Same as R
* Extract Col & Row
    * `a.col(1)` : Extract 1st column of `a` as Column matrix
    * `a.row(1)` : Extract 1st row of `a` as Row matrix

# Release 0.2.1 (2018-10-04)

* Update Documentation
    * Update README

# Release 0.2.0 (2018-10-04)

* Change structure
    * remove `ozone`
