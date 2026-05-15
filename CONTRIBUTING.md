# Contributing to Peroxide

Thank you for your interest in improving Peroxide!
Whether you want to report a bug, request a feature, or send a code change, this document explains the conventions the project follows.

## Reporting bugs and requesting features

Please open an issue at <https://github.com/Axect/Peroxide/issues>.

When reporting a bug, the following information makes it much easier to diagnose:

- A minimal Rust snippet that reproduces the problem
- Peroxide version (from `Cargo.toml` or `cargo tree | head`)
- Cargo features enabled (e.g. `O3-openblas`, `plot`, `nc`)
- Rust toolchain (`rustc --version`)
- Operating system and, when relevant, the BLAS / Python / HDF5 versions in use

Feature requests are also welcome through the issue tracker; describe the use case and, if you have one in mind, a suggested API shape.

## Asking for help

For questions that are not bug reports (usage, design discussion, "is this the right approach"), please use GitHub Discussions or open an issue tagged as a question.

## Submitting code changes

Peroxide follows the [Gitflow workflow]. A few practical rules:

1. **Do not commit to `master` directly.** `master` only receives merges from `dev` through reviewed pull requests.
2. **Branch from `dev`** for new work. Use the `features/<short-name>` convention (for example `features/dataframe-cleanup`). Open the pull request against `dev`.
3. **Run the test suite before pushing.** Enable the same Cargo features your change touches:
   ```sh
   cargo test --features "<features-you-touch>"
   ```
   For `O3` work, pick the convenience flag that matches your platform (`O3-openblas`, `O3-accelerate`, `O3-mkl`, or `O3-netlib`); see the README "Pre-requisite" section for details.
4. **Format and lint before opening the PR:**
   ```sh
   cargo fmt
   cargo clippy --all-targets
   ```
   The continuous integration on `push` / `pull_request` runs the default, `O3-openblas`, pure-Rust, and `plot` feature jobs, so keeping these green locally avoids round-trips.

## Source layout

- __src__
  - [lib.rs](src/lib.rs) : `mod` and `re-export`
  - __complex__: For complex vector, matrix & integrals.
    - [mod.rs](src/complex/mod.rs)
    - [integrate.rs](src/complex/integrate.rs) : Complex integral
    - [matrix.rs](src/complex/matrix.rs) : Complex matrix
    - [vector.rs](src/complex/vector.rs) : Complex vector
  - __fuga__ : Fuga for controlling numerical algorithms.
    - [mod.rs](src/fuga/mod.rs)
  - __macros__ : Macro files
    - [julia_macro.rs](src/macros/julia_macro.rs) : Julia like macro
    - [matlab_macro.rs](src/macros/matlab_macro.rs) : MATLAB like macro
    - [mod.rs](src/macros/mod.rs)
    - [r_macro.rs](src/macros/r_macro.rs) : R like macro
  - __ml__ : For machine learning (_Beta_)
    - [mod.rs](src/ml/mod.rs)
    - [reg.rs](src/ml/reg.rs) : Regression tools
  - __numerical__ : To do numerical things
    - [mod.rs](src/numerical/mod.rs)
    - [eigen.rs](src/numerical/eigen.rs) : Eigenvalue, Eigenvector algorithm
    - [integral.rs](src/numerical/integral.rs) : Numerical integration
    - [interp.rs](src/numerical/interp.rs) : Interpolation
    - [newton.rs](src/numerical/newton.rs) : Newton's Method
    - [ode.rs](src/numerical/ode.rs) : Main ODE solver with various algorithms
    - [optimize.rs](src/numerical/optimize.rs) : Non-linear regression
    - [root.rs](src/numerical/root.rs) : Root finding
    - [spline.rs](src/numerical/spline.rs) : Cubic spline, Cubic Hermite spline & B-Spline
    - [utils.rs](src/numerical/utils.rs) : Utils to do numerical things (e.g. jacobian)
  - __prelude__ : Prelude for using simple
    - [mod.rs](src/prelude/mod.rs)
    - [simpler.rs](src/prelude/simpler.rs) : Provides more simple api
  - __special__ : Special functions written in pure Rust (Wrapper of `puruspe`)
    - [mod.rs](src/special/mod.rs)
    - [function.rs](src/special/function.rs) : Special functions
  - __statistics__ : Statistical Tools
    - [mod.rs](src/statistics/mod.rs)
    - [dist.rs](src/statistics/dist.rs) : Probability distributions
    - [ops.rs](src/statistics/ops.rs) : Some probabilistic operations
    - [rand.rs](src/statistics/rand.rs) : Wrapper for `rand` crate & Piecewise Rejection Sampling
    - [stat.rs](src/statistics/stat.rs) : Statistical tools
  - __structure__ : Fundamental data structures
    - [mod.rs](src/structure/mod.rs)
    - [ad.rs](src/structure/ad.rs) : Automatic Differentiation (`Jet<N>` const-generic forward AD)
    - [dataframe.rs](src/structure/dataframe.rs) : Dataframe
    - [matrix.rs](src/structure/matrix.rs) : Matrix
    - [polynomial.rs](src/structure/polynomial.rs) : Polynomial
    - [sparse.rs](src/structure/sparse.rs) : For sparse structure (_Beta_)
    - [vector.rs](src/structure/vector.rs) : Extra tools for `Vec<f64>`
  - __traits__
    - [mod.rs](src/traits/mod.rs)
    - [fp.rs](src/traits/fp.rs) : Functional programming toolbox
    - [general.rs](src/traits/general.rs) : General algorithms
    - [math.rs](src/traits/math.rs) : Mathematics
    - [matrix.rs](src/traits/matrix.rs) : Matrix traits
    - [mutable.rs](src/traits/mutable.rs) : Mutable toolbox
    - [num.rs](src/traits/num.rs) : Number, Real and more operations
    - [pointer.rs](src/traits/pointer.rs) : Matrix pointer and Vector pointer for convenience
    - [stable.rs](src/traits/stable.rs) : Implement nightly-only features in stable
    - [sugar.rs](src/traits/sugar.rs) : Syntactic sugar for Vector
  - __util__
    - [mod.rs](src/util/mod.rs)
    - [api.rs](src/util/api.rs) : Matrix constructor for various language style
    - [low_level.rs](src/util/low_level.rs) : Low-level tools
    - [non_macro.rs](src/util/non_macro.rs) : Primordial version of macros
    - [plot.rs](src/util/plot.rs) : To draw plot (using `pyo3`)
    - [print.rs](src/util/print.rs) : To print conveniently
    - [useful.rs](src/util/useful.rs) : Useful utils to implement library
    - [wrapper.rs](src/util/wrapper.rs) : Wrapper for other crates (e.g. rand)
    - [writer.rs](src/util/writer.rs) : More convenient write system

## Code of conduct

By contributing you agree to abide by the project's
[Code of Conduct](CODE_OF_CONDUCT.md).

Thanks for all contributions!

[Gitflow workflow]: https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow
