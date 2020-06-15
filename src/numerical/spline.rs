use std::cmp::{max, min};
use std::ops::{Index, Range};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::convert::{From, Into};
#[allow(unused_imports)]
use crate::structure::matrix::*;
#[allow(unused_imports)]
use crate::structure::polynomial::*;
#[allow(unused_imports)]
use crate::structure::vector::*;
#[allow(unused_imports)]
use crate::util::non_macro::*;
use crate::traits::{
    num::PowOps,
};

/// Cubic Spline (Natural)
///
/// # Description
///
/// Implement traits of Natural cubic splines, Arne Morten Kvarving.
///
/// # Type
/// (Vec<f64>, Vec<f64>) -> Vec<Polynomial>
///
/// # Examples
/// ```
/// #[macro_use]
/// extern crate peroxide;
/// use peroxide::fuga::*;
///
/// fn main() {
///     let x = c!(0.9, 1.3, 1.9, 2.1);
///     let y = c!(1.3, 1.5, 1.85, 2.1);
///
///     let s = cubic_spline(x, y);
///
///     for i in 0 .. s.len() {
///         s[i].print();
///     }
///
///     // -0.2347x^3 + 0.6338x^2 - 0.0329x + 0.9873
///     // 0.9096x^3 - 3.8292x^2 + 5.7691x - 1.5268
///     // -2.2594x^3 + 14.2342x^2 - 28.5513x + 20.2094
/// }
/// ```
pub fn cubic_spline(node_x: Vec<f64>, node_y: Vec<f64>) -> Vec<Polynomial> {
    let spline = CubicSpline::from_nodes(node_x, node_y);
    spline.into()
}

/// Cubic Spline (Natural)
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CubicSpline {
    polynomials: Vec<(Range<f64>, Polynomial)>,
}

impl CubicSpline {
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let x = c!(0.9, 1.3, 1.9, 2.1);
    ///     let y = c!(1.3, 1.5, 1.85, 2.1);
    ///
    ///     let s = CubicSpline::from_nodes(x, y);
    ///
    ///     for i in 0 .. 4 {
    ///         println!("{}", s.eval(i as f64 / 2.0));
    ///     }
    /// }
    /// ```
    pub fn from_nodes(node_x: Vec<f64>, node_y: Vec<f64>) -> Self {
        let polynomials = CubicSpline::cubic_spline(&node_x, &node_y);
        CubicSpline {
            polynomials: CubicSpline::ranged(&node_x, polynomials),
        }
    }

    fn ranged(node_x: &Vec<f64>, polynomials: Vec<Polynomial>) -> Vec<(Range<f64>, Polynomial)> {
        let len = node_x.len();
        polynomials
            .into_iter()
            .enumerate()
            .filter(|(i, _)| i + 1 < len)
            .map(|(i, polynomial)| {
                (
                    Range {
                        start: node_x[i],
                        end: node_x[i + 1],
                    },
                    polynomial,
                )
            })
            .collect()
    }

    fn cubic_spline(node_x: &Vec<f64>, node_y: &Vec<f64>) -> Vec<Polynomial> {
        //! Pre calculated variables
        //! node_x: n+1
        //! node_y: n+1
        //! h     : n
        //! b     : n
        //! v     : n
        //! u     : n
        //! z     : n+1
        let n = node_x.len() - 1;
        assert_eq!(n, node_y.len() - 1);

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
        return s;
    }

    /// Evaluate cubic spline with value
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let x = c!(0.9, 1.3, 1.9, 2.1);
    ///     let y = c!(1.3, 1.5, 1.85, 2.1);
    ///
    ///     let s = CubicSpline::from_nodes(x, y);
    ///
    ///     s.eval(2.0);
    /// }
    /// ```
    pub fn eval<T>(&self, x: T) -> f64
    where
        T: std::convert::Into<f64> + Copy,
    {
        let x = x.into();

        self.polynomial(x).eval(x)
    }

    /// Returns a reference the `Polynomial` at the given point `x`.
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    ///
    /// fn main() {
    ///     let x = c!(0.9, 1.3, 1.9, 2.1);
    ///     let y = c!(1.3, 1.5, 1.85, 2.1);
    ///
    ///     let s = CubicSpline::from_nodes(x, y);
    ///
    ///     let p = s.polynomial(2.0);
    ///     let v = p.eval(1.9);
    ///
    ///     assert_eq!((v * 100.0).round() / 100.0, 1.85)
    /// }
    /// ```
    ///
    /// If `x` is outside of the range of polynomials, the first or last polynomial will be
    /// returned, depending if `x` is lower of the first interpolation point or higher of the last
    /// interpolation point.
    pub fn polynomial<T>(&self, x: T) -> &Polynomial
        where
        T: std::convert::Into<f64> + Copy,
    {
        let x  = x.into();

        let index = match self.polynomials.binary_search_by(|(range, _)| {
            if range.contains(&x) {
                core::cmp::Ordering::Equal
            } else if x < range.start {
                core::cmp::Ordering::Greater
            } else {
                core::cmp::Ordering::Less
            }
        }) {
            Ok(index) => index,
            Err(index) => max(0, min(index, self.polynomials.len() - 1)),
        };

        &self.polynomials[index].1
    }

    pub fn extend_with_nodes(&mut self, node_x: Vec<f64>, node_y: Vec<f64>) {
        let mut ext_node_x = Vec::with_capacity(node_x.len() + 1);
        let mut ext_node_y = Vec::with_capacity(node_x.len() + 1);

        let (r, polynomial) = &self.polynomials[self.polynomials.len() - 1];
        ext_node_x.push(r.end);
        ext_node_y.push(polynomial.eval(r.end));

        ext_node_x.extend(node_x);
        ext_node_y.extend(node_y);

        let polynomials = CubicSpline::ranged(
            &ext_node_x,
            CubicSpline::cubic_spline(&ext_node_x, &ext_node_y),
        );

        self.polynomials.extend(polynomials);
    }

    /// Returns the number of polynimials that describe the `CubicSpline`
    ///
    /// # Examples
    /// ```
    /// #[macro_use]
    /// extern crate peroxide;
    /// use peroxide::fuga::*;
    /// 
    /// fn main() {
    ///     let x = c!(0.9, 1.3, 1.9, 2.1);
    ///     let y = c!(1.3, 1.5, 1.85, 2.1);
    ///
    ///     let s = CubicSpline::from_nodes(x, y);
    ///
    ///     assert_eq!(s.number_of_polynomials(), 3);
    /// }
    /// ```
    pub fn number_of_polynomials(&self) -> usize {
        self.polynomials.len()
    }
}

impl std::convert::Into<Vec<Polynomial>> for CubicSpline {
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
    fn diff(&self) -> Self {
        let mut polynomials: Vec<(Range<f64>, Polynomial)> = self.clone().into();

        polynomials = polynomials
            .into_iter()
            .map(|(r, poly)| (r, poly.diff()))
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
}
