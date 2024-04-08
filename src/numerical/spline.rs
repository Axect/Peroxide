//! Spline interpolations
//!
//! # Available splines
//!
//! * Cubic spline
//! * Cubic Hermite spline
//!
//! # `Spline` trait
//!
//! ## Methods
//!
//! Let `T: Into<f64> + Copy`
//! * `fn eval<T>(&self, x: T) -> f64` : Evaluate the spline at x
//! * `fn eval_vec<T>(&self, v: &[T]) -> Vec<f64>` : Evaluate spline values for an array v
//! * `fn polynomial_at<T>(&self, x: T) -> &Polynomial` : Get the polynomial at x
//! * `fn number_of_polynomials(&self) -> usize` : Get the number of polynomials
//! * `fn get_ranged_polynomials(&self) -> &Vec<(Range<f64>, Polynomial)>` : Get the polynomials
//! * `fn eval_with_cond<F>(&self, x: f64, cond: F) -> f64` : Evaluate the spline at x, with a condition
//! * `fn eval_vec_with_cond<F>(&self, v: &[f64], cond: F) -> Vec<f64>` : Evaluate spline values for an array v, with a condition
//!
//! # Low-level interface
//!
//! ## Members
//!
//! * `CubicSpline`: Structure for cubic spline
//!     * `fn from_nodes(node_x: &[f64], node_y: &[f64]) -> Self` : Create a cubic spline from nodes
//!     * `fn extend_with_nodes(&mut self, node_x: Vec<f64>, node_y: Vec<f64>)` : Extend the spline with nodes
//! * `CubicHermiteSpline`: Structure for cubic Hermite spline
//!     * `fn from_nodes_with_slopes(node_x: &[f64], node_y: &[f64], m: &[f64]) -> Self` : Create a Cubic Hermite spline from nodes with slopes
//!     * `fn from_nodes(node_x: &[f64], node_y: &[f64], slope_method: SlopeMethod) -> Self` : Create a Cubic Hermite spline from nodes with slope estimation methods
//! * `SlopeMethod`: Enum for slope estimation methods
//!     * `Akima`: Akima's method to estimate slopes ([Akima (1970)](https://dl.acm.org/doi/abs/10.1145/321607.321609))
//!     * `Quadratic`: Using quadratic interpolation to estimate slopes
//!
//! ## Usage
//!
//! ```rust
//! use peroxide::fuga::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let x = seq(0, 10, 1);
//!     let y = x.fmap(|t| t.sin());
//!     
//!     let cs = CubicSpline::from_nodes(&x, &y)?;
//!     let cs_akima = CubicHermiteSpline::from_nodes(&x, &y, Akima)?;
//!     let cs_quad = CubicHermiteSpline::from_nodes(&x, &y, Quadratic)?;
//!
//!     cs.polynomial_at(0f64).print();
//!     cs_akima.polynomial_at(0f64).print();
//!     cs_quad.polynomial_at(0f64).print();
//!     // -0.1523x^3 + 0.9937x
//!     // 0.1259x^3 - 0.5127x^2 + 1.2283x
//!     // -0.0000x^3 - 0.3868x^2 + 1.2283x
//!
//!     let new_x = seq(4, 6, 0.1);
//!     let new_y = new_x.fmap(|t| t.sin());
//!
//!     let y_cs = cs.eval_vec(&new_x);
//!     let y_akima = cs_akima.eval_vec(&new_x);
//!     let y_quad = cs_quad.eval_vec(&new_x);
//!
//!     let mut df = DataFrame::new(vec![]);
//!     df.push("x", Series::new(new_x));
//!     df.push("y", Series::new(new_y));
//!     df.push("y_cs", Series::new(y_cs));
//!     df.push("y_akima", Series::new(y_akima));
//!     df.push("y_quad", Series::new(y_quad));
//!
//!     df.print();
//!     //          x       y    y_cs y_akima  y_quad
//!     //  r[0]    5 -0.9589 -0.9589 -0.9589 -0.9589
//!     //  r[1]  5.2 -0.8835 -0.8826 -0.8583 -0.8836
//!     //  r[2]  5.4 -0.7728 -0.7706 -0.7360 -0.7629
//!     //  r[3]  5.6 -0.6313 -0.6288 -0.5960 -0.6120
//!     //  r[4]  5.8 -0.4646 -0.4631 -0.4424 -0.4459
//!     //  r[5]    6 -0.2794 -0.2794 -0.2794 -0.2794
//!
//!     Ok(())
//! }
//! ```
//!
//! # High-level interface
//!
//! ## Functions
//!
//! * `fn cubic_spline(node_x: &[f64], node_y: &[f64]) -> CubicSpline` : Create a cubic spline from nodes
//! * `fn cubic_hermite_spline(node_x: &[f64], node_y: &[f64], m: &[f64]) -> CubicHermiteSpline` : Create a cubic Hermite spline from nodes with slopes
//!
//! ## Usage
//!
//! ```rust
//! use peroxide::fuga::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let x = seq(0, 10, 1);
//!     let y = x.fmap(|t| t.sin());
//!     
//!     let cs = cubic_spline(&x, &y)?;
//!     let cs_akima = cubic_hermite_spline(&x, &y, Akima)?;
//!     let cs_quad = cubic_hermite_spline(&x, &y, Quadratic)?;
//!
//!     cs.polynomial_at(0f64).print();
//!     cs_akima.polynomial_at(0f64).print();
//!     cs_quad.polynomial_at(0f64).print();
//!     // -0.1523x^3 + 0.9937x
//!     // 0.1259x^3 - 0.5127x^2 + 1.2283x
//!     // -0.0000x^3 - 0.3868x^2 + 1.2283x
//!
//!     let new_x = seq(4, 6, 0.1);
//!     let new_y = new_x.fmap(|t| t.sin());
//!
//!     let y_cs = cs.eval_vec(&new_x);
//!     let y_akima = cs_akima.eval_vec(&new_x);
//!     let y_quad = cs_quad.eval_vec(&new_x);
//!
//!     let mut df = DataFrame::new(vec![]);
//!     df.push("x", Series::new(new_x));
//!     df.push("y", Series::new(new_y));
//!     df.push("y_cs", Series::new(y_cs));
//!     df.push("y_akima", Series::new(y_akima));
//!     df.push("y_quad", Series::new(y_quad));
//!
//!     df.print();
//!     //          x       y    y_cs y_akima  y_quad
//!     //  r[0]    5 -0.9589 -0.9589 -0.9589 -0.9589
//!     //  r[1]  5.2 -0.8835 -0.8826 -0.8583 -0.8836
//!     //  r[2]  5.4 -0.7728 -0.7706 -0.7360 -0.7629
//!     //  r[3]  5.6 -0.6313 -0.6288 -0.5960 -0.6120
//!     //  r[4]  5.8 -0.4646 -0.4631 -0.4424 -0.4459
//!     //  r[5]    6 -0.2794 -0.2794 -0.2794 -0.2794
//!
//!     Ok(())
//! }
//! ```
//!
//! # Calculus with splines
//!
//! ## Usage
//!
//! ```rust
//! use peroxide::fuga::*;
//! use std::f64::consts::PI;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let x = seq(0, 10, 1);
//!     let y = x.fmap(|t| t.sin());
//!     
//!     let cs = cubic_spline(&x, &y)?;
//!     let cs_akima = cubic_hermite_spline(&x, &y, Akima)?;
//!     let cs_quad = cubic_hermite_spline(&x, &y, Quadratic)?;
//!
//!     println!("============ Polynomial at x=0 ============");
//!
//!     cs.polynomial_at(0f64).print();
//!     cs_akima.polynomial_at(0f64).print();
//!     cs_quad.polynomial_at(0f64).print();
//!
//!     println!("============ Derivative at x=0 ============");
//!
//!     cs.derivative().polynomial_at(0f64).print();
//!     cs_akima.derivative().polynomial_at(0f64).print();
//!     cs_quad.derivative().polynomial_at(0f64).print();
//!
//!     println!("============ Integral at x=0 ============");
//!
//!     cs.integral().polynomial_at(0f64).print();
//!     cs_akima.integral().polynomial_at(0f64).print();
//!     cs_quad.integral().polynomial_at(0f64).print();
//!
//!     println!("============ Integrate from x=0 to x=pi ============");
//!
//!     cs.integrate((0f64, PI)).print();
//!     cs_akima.integrate((0f64, PI)).print();
//!     cs_quad.integrate((0f64, PI)).print();
//!
//!     // ============ Polynomial at x=0 ============
//!     // -0.1523x^3 + 0.9937x
//!     // 0.1259x^3 - 0.5127x^2 + 1.2283x
//!     // -0.0000x^3 - 0.3868x^2 + 1.2283x
//!     // ============ Derivative at x=0 ============
//!     // -0.4568x^2 + 0.9937
//!     // 0.3776x^2 - 1.0254x + 1.2283
//!     // -0.0000x^2 - 0.7736x + 1.2283
//!     // ============ Integral at x=0 ============
//!     // -0.0381x^4 + 0.4969x^2
//!     // 0.0315x^4 - 0.1709x^3 + 0.6141x^2
//!     // -0.0000x^4 - 0.1289x^3 + 0.6141x^2
//!     // ============ Integrate from x=0 to x=pi ============
//!     // 1.9961861265456702
//!     // 2.0049920614062775
//!     // 2.004327391790717
//!
//!     Ok(())
//! }
//! ```
//!
//! # References
//!
//! * Gary D. Knott, *Interpolating Splines*, Birkhäuser Boston, MA, (2000).

use self::SplineError::{NotEnoughNodes, NotEqualNodes, NotEqualSlopes, RedundantNodeX};
#[allow(unused_imports)]
use crate::structure::matrix::*;
#[allow(unused_imports)]
use crate::structure::polynomial::*;
#[allow(unused_imports)]
use crate::structure::vector::*;
#[allow(unused_imports)]
use crate::util::non_macro::*;
use crate::util::useful::zip_range;
use peroxide_num::PowOps;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::cmp::{max, min};
use std::convert::{From, Into};
use std::ops::{Index, Range};
use thiserror::Error;

/// Trait for spline interpolation
///
/// # Available Splines
///
/// - `CubicSpline`
/// - `CubicHermiteSpline`
pub trait Spline {
    fn eval<T: Into<f64> + Copy>(&self, x: T) -> f64 {
        let x = x.into();

        self.polynomial_at(x).eval(x)
    }

    fn eval_vec<T: Into<f64> + Copy>(&self, v: &[T]) -> Vec<f64> {
        let mut result = vec![0f64; v.len()];

        for (i, x) in v.iter().enumerate() {
            result[i] = self.eval(*x);
        }

        result
    }

    fn polynomial_at<T: Into<f64> + Copy>(&self, x: T) -> &Polynomial {
        let x = x.into();

        let poly = self.get_ranged_polynomials();

        let index = match poly.binary_search_by(|(range, _)| {
            if range.contains(&x) {
                core::cmp::Ordering::Equal
            } else if x < range.start {
                core::cmp::Ordering::Greater
            } else {
                core::cmp::Ordering::Less
            }
        }) {
            Ok(index) => index,
            Err(index) => max(0, min(index, poly.len() - 1)),
        };

        &poly[index].1
    }

    fn number_of_polynomials(&self) -> usize {
        self.get_ranged_polynomials().len()
    }

    fn get_ranged_polynomials(&self) -> &Vec<(Range<f64>, Polynomial)>;

    fn eval_with_cond<F: Fn(f64) -> f64>(&self, x: f64, cond: F) -> f64 {
        cond(self.eval(x))
    }

    fn eval_vec_with_cond<F: Fn(f64) -> f64 + Copy>(&self, x: &[f64], cond: F) -> Vec<f64> {
        x.iter().map(|&x| self.eval_with_cond(x, cond)).collect()
    }
}

// =============================================================================
// High level functions
// =============================================================================
/// Cubic Spline (Natural)
///
/// # Description
///
/// Implement traits of Natural cubic splines, by Arne Morten Kvarving.
///
/// # Type
/// `(&[f64], &[f64]) -> Cubic Spline`
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let x = c!(0.9, 1.3, 1.9, 2.1);
///     let y = c!(1.3, 1.5, 1.85, 2.1);
///
///     let s = cubic_spline(&x, &y)?;
///
///     let new_x = c!(1, 1.5, 2.0);
///
///     // Generate Cubic polynomial
///     for t in new_x.iter() {
///         s.polynomial_at(*t).print();
///     }
///     // -0.2347x^3 + 0.6338x^2 - 0.0329x + 0.9873
///     // 0.9096x^3 - 3.8292x^2 + 5.7691x - 1.5268
///     // -2.2594x^3 + 14.2342x^2 - 28.5513x + 20.2094
///
///     // Evaluation
///     for t in new_x.iter() {
///         s.eval(*t).print();
///     }
///
///     Ok(())
/// }
/// ```
pub fn cubic_spline(node_x: &[f64], node_y: &[f64]) -> Result<CubicSpline, SplineError> {
    CubicSpline::from_nodes(node_x, node_y)
}

pub fn cubic_hermite_spline(
    node_x: &[f64],
    node_y: &[f64],
    slope_method: SlopeMethod,
) -> Result<CubicHermiteSpline, SplineError> {
    CubicHermiteSpline::from_nodes(node_x, node_y, slope_method)
}

// =============================================================================
// Cubic Spline
// =============================================================================
/// Cubic Spline (Natural)
///
/// # Description
///
/// Implement traits of Natural cubic splines, by Arne Morten Kvarving.
///
/// # Type
/// `(&[f64], &[f64]) -> Cubic Spline`
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let x = c!(0.9, 1.3, 1.9, 2.1);
///     let y = c!(1.3, 1.5, 1.85, 2.1);
///
///     let s = cubic_spline(&x, &y)?;
///
///     let new_x = c!(1, 1.5, 2.0);
///
///     // Generate Cubic polynomial
///     for t in new_x.iter() {
///         s.polynomial_at(*t).print();
///     }
///     // -0.2347x^3 + 0.6338x^2 - 0.0329x + 0.9873
///     // 0.9096x^3 - 3.8292x^2 + 5.7691x - 1.5268
///     // -2.2594x^3 + 14.2342x^2 - 28.5513x + 20.2094
///
///     // Evaluation
///     for t in new_x.iter() {
///         s.eval(*t).print();
///     }
///
///     Ok(())
/// }
/// ```
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CubicSpline {
    polynomials: Vec<(Range<f64>, Polynomial)>,
}

impl Spline for CubicSpline {
    fn get_ranged_polynomials(&self) -> &Vec<(Range<f64>, Polynomial)> {
        &self.polynomials
    }
}

#[derive(Debug, Copy, Clone, Error)]
pub enum SplineError {
    #[error("node_x has less than 3 elements")]
    NotEnoughNodes,
    #[error("node_x and node_y have different lengths")]
    NotEqualNodes,
    #[error("nodes and slopes have different lengths")]
    NotEqualSlopes,
    #[error("there are redundant nodes in node_x")]
    RedundantNodeX,
}

impl CubicSpline {
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let x = c!(0.9, 1.3, 1.9, 2.1);
    ///     let y = c!(1.3, 1.5, 1.85, 2.1);
    ///
    ///     let s = CubicSpline::from_nodes(&x, &y)?;
    ///
    ///     for i in 0 .. 4 {
    ///         println!("{}", s.eval(i as f64 / 2.0));
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn from_nodes(node_x: &[f64], node_y: &[f64]) -> Result<Self, SplineError> {
        let polynomials = CubicSpline::cubic_spline(node_x, node_y)?;
        Ok(CubicSpline {
            polynomials: zip_range(node_x, &polynomials),
        })
    }

    fn cubic_spline(node_x: &[f64], node_y: &[f64]) -> Result<Vec<Polynomial>, SplineError> {
        //! Pre calculated variables
        //! node_x: n+1
        //! node_y: n+1
        //! h     : n
        //! b     : n
        //! v     : n
        //! u     : n
        //! z     : n+1
        let n = node_x.len() - 1;
        if n < 2 {
            return Err(NotEnoughNodes);
        }
        if n != node_y.len() - 1 {
            return Err(NotEqualNodes);
        }

        // Pre-calculations
        let mut h = vec![0f64; n];
        let mut b = vec![0f64; n];
        let mut v = vec![0f64; n];
        let mut u = vec![0f64; n];
        for i in 0..n {
            if i == 0 {
                h[i] = node_x[i + 1] - node_x[i];
                b[i] = (node_y[i + 1] - node_y[i]) / h[i];
            } else {
                h[i] = node_x[i + 1] - node_x[i];
                b[i] = (node_y[i + 1] - node_y[i]) / h[i];
                v[i] = 2. * (h[i] + h[i - 1]);
                u[i] = 6. * (b[i] - b[i - 1]);
            }
        }

        // Tri-diagonal matrix
        let mut m = matrix(vec![0f64; (n - 1) * (n - 1)], n - 1, n - 1, Col);
        for i in 0..n - 2 {
            m[(i, i)] = v[i + 1];
            m[(i + 1, i)] = h[i + 1];
            m[(i, i + 1)] = h[i + 1];
        }
        m[(n - 2, n - 2)] = v[n - 1];

        // Calculate z
        let z_inner = m.inv() * Vec::from(&u[1..]);
        let mut z = vec![0f64];
        z.extend(&z_inner);
        z.push(0f64);

        // Declare empty spline
        let mut s: Vec<Polynomial> = Vec::new();

        // Main process
        for i in 0..n {
            // Memoization
            let t_i = node_x[i];
            let t_i1 = node_x[i + 1];
            let z_i = z[i];
            let z_i1 = z[i + 1];
            let h_i = h[i];
            let y_i = node_y[i];
            let y_i1 = node_y[i + 1];
            let temp1 = poly(vec![1f64, -t_i]);
            let temp2 = poly(vec![1f64, -t_i1]);

            let term1 = temp1.powi(3) * (z_i1 / (6f64 * h_i));
            let term2 = temp2.powi(3) * (-z_i / (6f64 * h_i));
            let term3 = temp1 * (y_i1 / h_i - z_i1 * h_i / 6.);
            let term4 = temp2 * (-y_i / h_i + h_i * z_i / 6.0);

            s.push(term1 + term2 + term3 + term4);
        }
        Ok(s)
    }

    /// Extends the spline with the given nodes.
    ///
    /// # Description
    ///
    /// The method ensures that the transition between each polynomial is smooth and that the spline
    /// interpolation of the new nodes is calculated around `x = 0` in order to avoid that
    /// successive spline extensions with large x values become inaccurate.
    pub fn extend_with_nodes(&mut self, node_x: Vec<f64>, node_y: Vec<f64>) -> Result<(), SplineError> {
        let mut ext_node_x = Vec::with_capacity(node_x.len() + 1);
        let mut ext_node_y = Vec::with_capacity(node_x.len() + 1);

        let (r, polynomial) = &self.polynomials[self.polynomials.len() - 1];
        ext_node_x.push(0.0f64);
        ext_node_y.push(polynomial.eval(r.end));

        // translating the node towards x = 0 increases accuracy tremendously
        let tx = r.end;
        ext_node_x.extend(node_x.into_iter().map(|x| x - tx));
        ext_node_y.extend(node_y);

        let polynomials = zip_range(
            &ext_node_x,
            &CubicSpline::cubic_spline(&ext_node_x, &ext_node_y)?,
        );

        self.polynomials
            .extend(polynomials.into_iter().map(|(r, p)| {
                (
                    Range {
                        start: r.start + tx,
                        end: r.end + tx,
                    },
                    p.translate_x(tx),
                )
            }));

        Ok(())
    }
}

impl Into<Vec<Polynomial>> for CubicSpline {
    fn into(self) -> Vec<Polynomial> {
        self.polynomials
            .into_iter()
            .map(|(_, polynomial)| polynomial)
            .collect()
    }
}

impl From<Vec<(Range<f64>, Polynomial)>> for CubicSpline {
    fn from(polynomials: Vec<(Range<f64>, Polynomial)>) -> Self {
        CubicSpline { polynomials }
    }
}

impl Into<Vec<(Range<f64>, Polynomial)>> for CubicSpline {
    fn into(self) -> Vec<(Range<f64>, Polynomial)> {
        self.polynomials
    }
}

impl Index<usize> for CubicSpline {
    type Output = (Range<f64>, Polynomial);

    fn index(&self, index: usize) -> &Self::Output {
        &self.polynomials[index]
    }
}

impl Calculus for CubicSpline {
    fn derivative(&self) -> Self {
        let mut polynomials: Vec<(Range<f64>, Polynomial)> = self.clone().into();

        polynomials = polynomials
            .into_iter()
            .map(|(r, poly)| (r, poly.derivative()))
            .collect();

        Self::from(polynomials)
    }

    fn integral(&self) -> Self {
        let mut polynomials: Vec<(Range<f64>, Polynomial)> = self.clone().into();

        polynomials = polynomials
            .into_iter()
            .map(|(r, poly)| (r, poly.integral()))
            .collect();

        Self::from(polynomials)
    }

    fn integrate<T: Into<f64> + Copy>(&self, interval: (T, T)) -> f64 {
        let (a, b) = interval;
        let a = a.into();
        let b = b.into();

        let mut s = 0f64;
        for (r, p) in self.polynomials.iter() {
            if r.start > b {
                break;
            } else if r.end < a {
                continue;
            } else {
                // r.start <= b, r.end >= a
                let x = r.start.max(a);
                let y = r.end.min(b);
                s += p.integrate((x, y));
            }
        }
        s
    }
}

// =============================================================================
// Cubic Hermite Spline
// =============================================================================
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CubicHermiteSpline {
    polynomials: Vec<(Range<f64>, Polynomial)>,
}

impl Spline for CubicHermiteSpline {
    fn get_ranged_polynomials(&self) -> &Vec<(Range<f64>, Polynomial)> {
        &self.polynomials
    }
}

impl CubicHermiteSpline {
    pub fn from_nodes_with_slopes(
        node_x: &[f64],
        node_y: &[f64],
        m: &[f64],
    ) -> Result<Self, SplineError> {
        let n = node_x.len();
        if n < 3 {
            return Err(NotEnoughNodes);
        }
        if n != node_y.len() {
            return Err(NotEqualNodes);
        }
        if n != m.len() {
            return Err(NotEqualSlopes);
        }

        let mut r = vec![Range::default(); node_x.len() - 1];
        let mut u = vec![Polynomial::default(); node_x.len() - 1];

        for i in 0..node_x.len() - 1 {
            let a_i = node_y[i];
            let b_i = m[i];
            let dx = node_x[i + 1] - node_x[i];
            let dy = node_y[i + 1] - node_y[i];
            let c_i = (3f64 * dy / dx - 2f64 * m[i] - m[i + 1]) / dx;
            let d_i = (m[i] + m[i + 1] - 2f64 * dy / dx) / dx.powi(2);

            let p = Polynomial::new(vec![1f64, -node_x[i]]);

            r[i] = Range {
                start: node_x[i],
                end: node_x[i + 1],
            };
            u[i] = p.powi(3) * d_i + p.powi(2) * c_i + p.clone() * b_i;
            u[i].coef[3] += a_i;
        }

        Ok(CubicHermiteSpline {
            polynomials: r.into_iter().zip(u).collect(),
        })
    }

    pub fn from_nodes(
        node_x: &[f64],
        node_y: &[f64],
        slope_method: SlopeMethod,
    ) -> Result<Self, SplineError> {
        match slope_method {
            SlopeMethod::Akima => CubicHermiteSpline::from_nodes_with_slopes(
                node_x,
                node_y,
                &akima_slopes(node_x, node_y)?,
            ),
            SlopeMethod::Quadratic => CubicHermiteSpline::from_nodes_with_slopes(
                node_x,
                node_y,
                &quadratic_slopes(node_x, node_y)?,
            ),
        }
    }
}

impl Into<Vec<Polynomial>> for CubicHermiteSpline {
    fn into(self) -> Vec<Polynomial> {
        self.polynomials
            .into_iter()
            .map(|(_, polynomial)| polynomial)
            .collect()
    }
}

impl From<Vec<(Range<f64>, Polynomial)>> for CubicHermiteSpline {
    fn from(polynomials: Vec<(Range<f64>, Polynomial)>) -> Self {
        CubicHermiteSpline { polynomials }
    }
}

impl Into<Vec<(Range<f64>, Polynomial)>> for CubicHermiteSpline {
    fn into(self) -> Vec<(Range<f64>, Polynomial)> {
        self.polynomials
    }
}

impl Index<usize> for CubicHermiteSpline {
    type Output = (Range<f64>, Polynomial);

    fn index(&self, index: usize) -> &Self::Output {
        &self.polynomials[index]
    }
}

impl Calculus for CubicHermiteSpline {
    fn derivative(&self) -> Self {
        let mut polynomials: Vec<(Range<f64>, Polynomial)> = self.clone().into();

        polynomials = polynomials
            .into_iter()
            .map(|(r, poly)| (r, poly.derivative()))
            .collect();

        Self::from(polynomials)
    }

    fn integral(&self) -> Self {
        let mut polynomials: Vec<(Range<f64>, Polynomial)> = self.clone().into();

        polynomials = polynomials
            .into_iter()
            .map(|(r, poly)| (r, poly.integral()))
            .collect();

        Self::from(polynomials)
    }

    fn integrate<T: Into<f64> + Copy>(&self, interval: (T, T)) -> f64 {
        let (a, b) = interval;
        let a = a.into();
        let b = b.into();

        let mut s = 0f64;
        for (r, p) in self.polynomials.iter() {
            if r.start > b {
                break;
            } else if r.end < a {
                continue;
            } else {
                // r.start <= b, r.end >= a
                let x = r.start.max(a);
                let y = r.end.min(b);
                s += p.integrate((x, y));
            }
        }
        s
    }
}

// =============================================================================
// Estimate Slopes
// =============================================================================
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SlopeMethod {
    Akima,
    Quadratic,
}

fn akima_slopes(x: &[f64], y: &[f64]) -> Result<Vec<f64>, SplineError> {
    if x.len() < 3 {
        return Err(NotEnoughNodes);
    }

    let mut m = vec![0f64; x.len()];
    let mut s = vec![0f64; x.len() + 3]; // -2, -1, 0, ..., x.len()-1, x.len()

    let l_i = lagrange_polynomial(x[0..3].to_vec(), y[0..3].to_vec());
    let l_f = lagrange_polynomial(x[x.len() - 3..].to_vec(), y[y.len() - 3..].to_vec());

    let x_i = x[0] - (x[1] - x[0]) / 10f64;
    let x_ii = x_i - (x[1] - x[0]) / 10f64;
    let x_f = x[x.len() - 1] + (x[x.len() - 1] - x[x.len() - 2]) / 10f64;
    let x_ff = x_f + (x[x.len() - 1] - x[x.len() - 2]) / 10f64;

    let y_i = l_i.eval(x_i);
    let y_ii = l_i.eval(x_ii);
    let y_f = l_f.eval(x_f);
    let y_ff = l_f.eval(x_ff);

    let new_x = concat(&concat(&vec![x_ii, x_i], &x.to_vec()), &vec![x_f, x_ff]);
    let new_y = concat(&concat(&vec![y_ii, y_i], &y.to_vec()), &vec![y_f, y_ff]);

    for i in 0..new_x.len() - 1 {
        let dx = new_x[i + 1] - new_x[i];
        if dx == 0f64 {
            return Err(RedundantNodeX);
        }
        s[i] = (new_y[i + 1] - new_y[i]) / dx;
    }

    for i in 0..x.len() {
        let j = i + 2;
        let ds_f = (s[j + 1] - s[j]).abs();
        let ds_i = (s[j - 1] - s[j - 2]).abs();

        m[i] = if ds_f == 0f64 && ds_i == 0f64 {
            (s[j - 1] + s[j]) / 2f64
        } else {
            (ds_f * s[j - 1] + ds_i * s[j]) / (ds_f + ds_i)
        };
    }
    Ok(m)
}

fn quadratic_slopes(x: &[f64], y: &[f64]) -> Result<Vec<f64>, SplineError> {
    if x.len() < 3 {
        return Err(NotEnoughNodes);
    }
    let mut m = vec![0f64; x.len()];
    let q_i = lagrange_polynomial(x[0..3].to_vec(), y[0..3].to_vec());
    let q_f = lagrange_polynomial(x[x.len() - 3..].to_vec(), y[y.len() - 3..].to_vec());

    m[0] = q_i.derivative().eval(x[0]);
    m[x.len() - 1] = q_f.derivative().eval(x[x.len() - 1]);

    for i in 1..x.len() - 1 {
        let q = lagrange_polynomial(x[i - 1..i + 2].to_vec(), y[i - 1..i + 2].to_vec());
        m[i] = q.derivative().eval(x[i]);
    }

    Ok(m)
}
