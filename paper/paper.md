---
title: "Peroxide: A Batteries-Included Numerical Computing Library for Rust"
tags:
  - Rust
  - numerical computing
  - automatic differentiation
  - ordinary differential equations
  - scientific computing
  - spline interpolation
authors:
  - name: Tae-Geun Kim
    orcid: 0009-0000-4229-2935
    affiliation: "1, 2"
affiliations:
  - name: Key Laboratory of Nuclear Physics and Ion-beam Application (MOE), Institute of Modern Physics, Fudan University, Shanghai 200433, China
    index: 1
  - name: RIKEN Center for Interdisciplinary Theoretical and Mathematical Sciences (iTHEMS), Wako, Saitama 351-0198, Japan
    index: 2
date: 21 March 2026
bibliography: paper.bib
---

# Summary

Peroxide is a numerical computing library for the Rust programming language that integrates eight domains within a single crate: linear algebra, ordinary differential equation (ODE) solving, adaptive quadrature, spline interpolation with exact polynomial calculus, root-finding, forward-mode automatic differentiation (AD), statistics, and multi-format DataFrame I/O.
We designed Peroxide for applied scientists, engineers, and students who need a cohesive numerical environment in Rust without assembling a multi-crate dependency graph.
A dual-module API provides progressive disclosure: the `prelude` module offers sensible defaults, while the `fuga` module exposes explicit algorithm selection, addressing Rust's lack of default function arguments.
Peroxide is released under the MIT/Apache-2.0 dual license and archived on Zenodo [@peroxide].

# Statement of need

Scientific computing workflows rarely confine themselves to a single numerical domain.
A typical pipeline may solve an ODE, interpolate the solution onto an irregular grid, compute a definite integral, find a root of a derived quantity, and serialize the results for statistical analysis.
In Python, SciPy [@scipy] and NumPy [@numpy] consolidate these capabilities under one namespace.
In Rust, they are scattered across specialized crates, each with its own type conventions and trait hierarchies.

\autoref{tab:comparison} surveys the Rust numerical ecosystem.
Specialist crates excel individually: diffsol [@diffsol] provides ODE/DAE solvers with sensitivity analysis; faer [@faer] offers pure-Rust linear algebra competitive with LAPACK; nalgebra [@nalgebra] provides matrix abstractions widely adopted across the ecosystem; argmin [@argmin] supplies over 20 optimization algorithms; Polars [@polars] dominates tabular data processing.
russell [@russell] covers the broadest range among alternatives (ODE solvers, quadrature, root-finding, and statistics) but lacks AD, spline calculus, and DataFrame support.
No single crate spans all eight domains shown in \autoref{tab:comparison}.
Peroxide fills this integration gap.

: Feature coverage of Rust numerical libraries. Each entry reflects the crate's own API, including optional feature flags; a dash indicates the feature is absent. Crates providing only N-dimensional array abstractions without numerical algorithms (e.g., ndarray) are omitted. \label{tab:comparison}

|                     | Peroxide | diffsol | russell | faer | nalgebra | argmin | num-dual | Polars |
|---------------------|:--------:|:-------:|:-------:|:----:|:--------:|:------:|:--------:|:------:|
| Linear algebra      | $\checkmark$        | —       | $\checkmark$       | $\checkmark$    | $\checkmark$        | —      | —        | —      |
| ODE solvers         | $\checkmark$        | $\checkmark$       | $\checkmark$       | —    | —        | —      | —        | —      |
| Quadrature          | $\checkmark$        | —       | $\checkmark$       | —    | —        | —      | —        | —      |
| Spline interp.      | $\checkmark$        | —       | —       | —    | —        | —      | —        | —      |
| Root-finding        | $\checkmark$        | —       | $\checkmark$       | —    | —        | $\checkmark$      | —        | —      |
| Automatic diff.     | $\checkmark$        | $\checkmark$       | —       | —    | —        | —      | $\checkmark$        | —      |
| Statistics          | $\checkmark$        | —       | $\checkmark$       | —    | —        | —      | —        | —      |
| DataFrame I/O       | $\checkmark$        | —       | —       | —    | —        | —      | —        | $\checkmark$      |

We do not claim superiority over any specialist crate within its domain.
Peroxide's contribution is the integration itself, enabled by design abstractions that connect domains: a `ButcherTableau` trait with blanket implementation that unifies all Runge-Kutta methods through compile-time constants, a `Calculus` trait that supports exact differentiation and integration of piecewise polynomial splines, const-generic dimensional typing in the root-finding interface, and a `Real` trait that allows AD-derived derivatives to flow into optimization and root-finding routines without glue code.

# State of the field

In automatic differentiation specifically, a review of available Rust AD crates on crates.io indicates that, to our knowledge, Peroxide is the only library combining const-generic derivative order, normalized Taylor coefficient storage, and true Taylor-mode propagation with $O(N^2)$ cost per elementary operation.
Most alternatives provide dual numbers up to second or third order with fixed type hierarchies [@numdual] or nested generics with exponential cost at higher orders [@autodiff_elrnv].
The ad-trait crate [@adtrait] offers both forward and reverse modes but is limited to first-order derivatives.

For ODE solving, diffsol [@diffsol] is the most feature-rich Rust crate, providing BDF and SDIRK methods with sensitivity analysis.
Peroxide's `ButcherTableau` trait-based architecture is a complementary approach: it encodes Runge-Kutta methods as compile-time constants, making it straightforward for users to add new methods by defining only coefficient arrays.

# Software design

**Architecture.**
We store matrices as a flat `Vec<f64>` with a `Shape` enum (`Row`/`Col`) that controls logical layout without copying data.
This design trades the rich type-level dimensionality of nalgebra for a memory model that maps directly to both the pure-Rust `matrixmultiply` crate [@matrixmultiply] and, when the `O3` feature flag is enabled, to OpenBLAS, passing raw pointers with stride and transpose flags in both cases without intermediate type conversions.
A `Real` trait abstracts over `f64` and `AD` (= `Jet<2>`), so the same function can compute both values and derivatives.

**Automatic differentiation.**
The const-generic type `Jet<N>` stores the function value $c_0 = f(a)$ and $N$ normalized Taylor coefficients $c_k = f^{(k)}(a)/k!$ [@griewank2008].
Multiplication follows the Cauchy product of truncated power series:

$$c_n(f \cdot g) = \sum_{i=0}^{n} c_i(f)\, c_{n-i}(g)$$

Because the coefficients are pre-divided by $k!$, no binomial coefficients appear, yielding $O(N^2)$ Taylor-mode propagation that avoids factorial overflow at higher orders.
Division and transcendental functions use analogous recurrence relations on normalized coefficients.
The `#[ad_function]` procedural macro (from the companion `peroxide-ad` crate) transforms a plain `fn(f64) -> f64` into versions accepting `Jet<N>`, automatically generating first- and second-derivative functions.

**ODE solvers.**
A `ButcherTableau` trait [@butcher2016] declares Runge-Kutta stage coefficients as compile-time constant arrays.
A blanket implementation `impl<BU: ButcherTableau> ODEIntegrator for BU` provides the stepping logic, so adding a new method requires only four constant slices.
We ship 10 integrators: four fixed-step explicit methods (RALS3, RK4, RALS4, RK5), five embedded adaptive pairs (BS23 [@bogacki1989], RKF45 [@fehlberg1969], DP45 [@dormand1980], TSIT45 [@tsitouras2011], RKF78 [@fehlberg1969]), and one implicit Gauss-Legendre method (GL4).

**Spline calculus.**
Cubic, cubic Hermite [@akima1970], and B-spline [@knott2000] interpolation are provided.
The `PolynomialSpline` trait exposes `polynomial_at(x)`, returning the `Polynomial` governing each interval.
The `Calculus` trait then operates symbolically: `derivative()` returns a new spline of analytic derivatives, and `integrate((a, b))` evaluates the exact definite integral by antidifferentiation, not finite differencing.
To our knowledge, no other Rust spline crate [@splines_crate; @bspline_crate] exposes piecewise polynomial objects for exact calculus.

**Root-finding and quadrature.**
`RootFindingProblem<const I, const O, T>` enforces input/output dimensions at compile time; the Jacobian `[[f64; I]; O]` is stack-allocated.
Five solvers are provided, including Broyden's method for multivariate systems.
Gauss-Kronrod adaptive quadrature [@laurie1997; @piessens1983] is available in six orders (G7K15 through G30K61), each with a relative-error variant.

**Statistics and I/O.**
Distributions (uniform, normal, gamma, beta, Student's $t$, and others) share a common trait with `sample`, `pdf`, and `cdf` methods.
The `DataFrame` structure unifies CSV, NetCDF [@netcdf], and Parquet [@parquet] serialization.

# Example

The progressive disclosure API allows the same integration task to be expressed at two levels of control:

```rust
// With prelude: sensible default (Gauss-Kronrod G7K15R, tol = 1e-4)
use peroxide::prelude::*;
let area = integrate(|x| x.sin(), (0.0, std::f64::consts::PI));
```

```rust
// With fuga: explicit algorithm selection
use peroxide::fuga::*;
use std::f64::consts::PI;
let area = integrate(|x| x.sin(), (0.0, PI), GaussLegendre(15));
```

# Research impact statement

Peroxide has been publicly available on crates.io since 2018 and has accumulated over 1,050,000 total downloads (as reported by crates.io) and 677 GitHub stars.
Independent researchers have adopted Peroxide in peer-reviewed work: @comitini2025 used Peroxide's Gauss-Kronrod quadrature as the numerical integration backend for their QCD screened massive expansion calculations, and @steuteville2024 identified Peroxide as one of the three most widely used Rust ODE solver packages in a survey conducted at the U.S. National Renewable Energy Laboratory.
The author's own research has also relied on Peroxide: @neural_hamilton employed its implicit Gauss-Legendre ODE solver for symplectic reference solutions, @hyperboliclr used it for dataset generation, and @jho2026 used its Gauss-Kronrod quadrature (G7K15R) and cubic Hermite spline interpolation for primordial black hole axion spectral computations.
The Zenodo archive [@peroxide] ensures long-term citability and reproducibility.

# Validation

We maintain a test suite of 26 test modules.
The AD module alone contains 115 unit tests verifying `Jet<N>` at orders $N = 1$ through $10$ against symbolic reference values.
Across all tested functions and derivative orders, `Jet<N>` achieves relative errors at the level of machine epsilon ($\sim 10^{-15}$); in contrast, central finite-difference approximations degrade to relative errors exceeding $10^{0}$ by order four due to the fundamental tension between truncation error $O(h^p)$ and cancellation error $O(\varepsilon_{\mathrm{mach}} / h^k)$ [@press2007].
ODE integrators are tested against problems with known analytic solutions.
Continuous integration runs via GitHub Actions on every push.
The repository provides 41 standalone examples, and API documentation with KaTeX-rendered mathematics is hosted at [https://axect.github.io/Peroxide_Doc](https://axect.github.io/Peroxide_Doc).
Community guidelines and contribution instructions are available in the repository's `CONTRIBUTING.md`.

# Limitations

We do not claim best-in-class performance in any single domain: diffsol provides richer DAE support, faer delivers faster dense factorizations, and argmin offers a broader optimization catalog.
The `Real` trait currently covers only `f64` and `Jet<2>`; extending it to arbitrary `Jet<N>` and generalizing the `ODEProblem` interface to generic types are planned.
We support forward-mode AD only; reverse-mode differentiation and interoperability with Rust's experimental `std::autodiff` are longer-term goals.
Sparse matrices, PDE solvers, and FFT lie outside the current scope.

# AI usage disclosure

A generative AI tool (Claude Opus 4, Anthropic, March 2026) was used for grammar checking and copy-editing of the manuscript text.
All scientific content, software design decisions, technical claims, and code were produced entirely by the human author.
The final manuscript was reviewed, edited, and approved by the author.

# Acknowledgements

We thank all contributors to the Peroxide repository and the Rust scientific computing community for their feedback.

# References
