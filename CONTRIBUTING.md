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

A high-level map of `src/`; see each module's `mod.rs` and the [API docs](https://docs.rs/peroxide) for details.

| Module | Purpose |
| ------------------ | ------------------------------------------------------------------ |
| [`structure`](src/structure) | Core data structures: `Matrix`, `Vec<f64>` extensions, `DataFrame`, `Polynomial`, `Jet<N>` forward AD |
| [`numerical`](src/numerical) | Numerical algorithms: ODE solvers, integration, interpolation, splines, root finding, optimization, eigenvalues |
| [`statistics`](src/statistics) | Probability distributions, RNG wrappers, ordered statistics |
| [`complex`](src/complex) | Complex vectors, matrices, and integrals (`complex` feature) |
| [`special`](src/special) | Special functions (wrapper of `puruspe`) |
| [`traits`](src/traits) | Shared trait definitions (math, functional programming, pointers) |
| [`macros`](src/macros) | R / MATLAB / Julia style macros |
| [`fuga`](src/fuga), [`prelude`](src/prelude) | The two user-facing import styles (explicit vs simple) |
| [`util`](src/util) | Constructors, printing, plotting, low-level helpers |
| [`ml`](src/ml) | Basic machine learning tools (beta) |
| [`grave`](src/grave) | Retired implementations kept for reference; not compiled |

## Code of conduct

By contributing you agree to abide by the project's
[Code of Conduct](CODE_OF_CONDUCT.md).

Thanks for all contributions!

[Gitflow workflow]: https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow
