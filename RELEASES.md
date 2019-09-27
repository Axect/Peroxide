# Release 0.17.0 (2019-09-27) (Candidates)

## DataFrame (Alpha stage)

### Added

* `structure/dataframe.rs`

### Implement level for DataFrame

- [x] Pretty print
- [x] Convert to Row matrix
- [ ] Convert to Col matrix
- [x] Call by header (`Index` with header)
- [ ] Call by row
- [ ] `IndexMut` with header
- [ ] Any column operations
- [ ] Read from csv
- [ ] Write to csv
- [ ] Plot

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
    use peroxide::*;

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
