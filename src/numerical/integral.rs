use crate::structure::polynomial::{lagrange_polynomial, Calculus, Polynomial};
use crate::traits::fp::FPVector;
use crate::util::non_macro::seq;

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Integral {
    GaussLegendre(usize),
    NewtonCotes(usize),
    G7K15(f64, u32),
    G10K21(f64, u32),
    G15K31(f64, u32),
    G20K41(f64, u32),
    G25K51(f64, u32),
    G30K61(f64, u32),
    G7K15R(f64, u32),
    G10K21R(f64, u32),
    G15K31R(f64, u32),
    G20K41R(f64, u32),
    G25K51R(f64, u32),
    G30K61R(f64, u32),
}

impl Integral {
    pub fn get_num_node(&self) -> usize {
        match self {
            Integral::GaussLegendre(n) => *n,
            Integral::NewtonCotes(n) => *n,
            _ => panic!("This method does not have a fixed number of nodes."),
        }
    }

    pub fn get_tol(&self) -> f64 {
        match self {
            Integral::G7K15(tol, _) => *tol,
            Integral::G10K21(tol, _) => *tol,
            Integral::G15K31(tol, _) => *tol,
            Integral::G20K41(tol, _) => *tol,
            Integral::G25K51(tol, _) => *tol,
            Integral::G30K61(tol, _) => *tol,
            Integral::G7K15R(tol, _) => *tol,
            Integral::G10K21R(tol, _) => *tol,
            Integral::G15K31R(tol, _) => *tol,
            Integral::G20K41R(tol, _) => *tol,
            Integral::G25K51R(tol, _) => *tol,
            Integral::G30K61R(tol, _) => *tol,
            _ => panic!("This method does not have a tolerance."),
        }
    }

    pub fn get_max_iter(&self) -> u32 {
        match self {
            Integral::G7K15(_, max_iter) => *max_iter,
            Integral::G10K21(_, max_iter) => *max_iter,
            Integral::G15K31(_, max_iter) => *max_iter,
            Integral::G20K41(_, max_iter) => *max_iter,
            Integral::G25K51(_, max_iter) => *max_iter,
            Integral::G30K61(_, max_iter) => *max_iter,
            Integral::G7K15R(_, max_iter) => *max_iter,
            Integral::G10K21R(_, max_iter) => *max_iter,
            Integral::G15K31R(_, max_iter) => *max_iter,
            Integral::G20K41R(_, max_iter) => *max_iter,
            Integral::G25K51R(_, max_iter) => *max_iter,
            Integral::G30K61R(_, max_iter) => *max_iter,
            _ => panic!("This method does not have a maximum number of iterations."),
        }
    }

    pub fn get_gauss_kronrod_order(&self) -> (u8, u8) {
        match self {
            Integral::G7K15(_, _) => (7, 15),
            Integral::G10K21(_, _) => (10, 21),
            Integral::G15K31(_, _) => (15, 31),
            Integral::G20K41(_, _) => (20, 41),
            Integral::G25K51(_, _) => (25, 51),
            Integral::G30K61(_, _) => (30, 61),
            Integral::G7K15R(_, _) => (7, 15),
            Integral::G10K21R(_, _) => (10, 21),
            Integral::G15K31R(_, _) => (15, 31),
            Integral::G20K41R(_, _) => (20, 41),
            Integral::G25K51R(_, _) => (25, 51),
            Integral::G30K61R(_, _) => (30, 61),
            _ => panic!("This method does not have a Gauss-Kronrod order."),
        }
    }

    pub fn is_relative(&self) -> bool {
        match self {
            Integral::G7K15R(_, _) => true,
            Integral::G10K21R(_, _) => true,
            Integral::G15K31R(_, _) => true,
            Integral::G20K41R(_, _) => true,
            Integral::G25K51R(_, _) => true,
            Integral::G30K61R(_, _) => true,
            _ => false,
        }
    }

    pub fn change_tol(&self, tol: f64) -> Self {
        match self {
            Integral::G7K15(_, max_iter) => Integral::G7K15(tol, *max_iter),
            Integral::G10K21(_, max_iter) => Integral::G10K21(tol, *max_iter),
            Integral::G15K31(_, max_iter) => Integral::G15K31(tol, *max_iter),
            Integral::G20K41(_, max_iter) => Integral::G20K41(tol, *max_iter),
            Integral::G25K51(_, max_iter) => Integral::G25K51(tol, *max_iter),
            Integral::G30K61(_, max_iter) => Integral::G30K61(tol, *max_iter),
            Integral::G7K15R(_, max_iter) => Integral::G7K15R(tol, *max_iter),
            Integral::G10K21R(_, max_iter) => Integral::G10K21R(tol, *max_iter),
            Integral::G15K31R(_, max_iter) => Integral::G15K31R(tol, *max_iter),
            Integral::G20K41R(_, max_iter) => Integral::G20K41R(tol, *max_iter),
            Integral::G25K51R(_, max_iter) => Integral::G25K51R(tol, *max_iter),
            Integral::G30K61R(_, max_iter) => Integral::G30K61R(tol, *max_iter),
            _ => panic!("This method does not have a tolerance."),
        }
    }

    pub fn change_max_iter(&self, max_iter: u32) -> Self {
        match self {
            Integral::G7K15(tol, _) => Integral::G7K15(*tol, max_iter),
            Integral::G10K21(tol, _) => Integral::G10K21(*tol, max_iter),
            Integral::G15K31(tol, _) => Integral::G15K31(*tol, max_iter),
            Integral::G20K41(tol, _) => Integral::G20K41(*tol, max_iter),
            Integral::G25K51(tol, _) => Integral::G25K51(*tol, max_iter),
            Integral::G30K61(tol, _) => Integral::G30K61(*tol, max_iter),
            Integral::G7K15R(tol, _) => Integral::G7K15R(*tol, max_iter),
            Integral::G10K21R(tol, _) => Integral::G10K21R(*tol, max_iter),
            Integral::G15K31R(tol, _) => Integral::G15K31R(*tol, max_iter),
            Integral::G20K41R(tol, _) => Integral::G20K41R(*tol, max_iter),
            Integral::G25K51R(tol, _) => Integral::G25K51R(*tol, max_iter),
            Integral::G30K61R(tol, _) => Integral::G30K61R(*tol, max_iter),
            _ => panic!("This method does not have a maximum number of iterations."),
        }
    }
}

/// Type that can be integrated using Newton Cotes Quadrature
pub trait NCIntegrable: std::ops::Sub<Self, Output = Self> + Sized {
    type NodeY;
    type NCPolynomial;

    /// Returns the image of `node_x` under function `f`, in a representation
    /// (`Self::NodeY`) suitable for separately computing one Lagrange polynomial
    /// for each of the degrees of freedom of `Self` (e.g. a `Vec` with as many
    /// entries as the number of degrees of freedom).
    fn compute_node_y<F>(f: F, node_x: &[f64]) -> Self::NodeY
    where
        F: Fn(f64) -> Self;

    /// Computes one [`Lagrange polynomial`](lagrange_polynomial) for each one
    /// of the degrees of freedom of `Self`, returning them in a representation
    /// (`Self::NCPolynomial`) suitable to separately performing operation on
    /// them (e.g. [`Vec<Polynomial>`](Polynomial)).
    fn compute_polynomial(node_x: &[f64], node_y: &Self::NodeY) -> Self::NCPolynomial;

    /// Separately integrates each of the polynomials obtained using
    /// [`compute_polynomial`](NCIntegrable::compute_polynomial).
    fn integrate_polynomial(p: &Self::NCPolynomial) -> Self::NCPolynomial;

    /// Separately evaluates each of the polynomial integrated using
    /// [`integrate_polynomial`](NCIntegrable::integrate_polynomial) at `x`,
    /// then recombines the result into a `Self`.
    fn evaluate_polynomial(p: &Self::NCPolynomial, x: f64) -> Self;
}

impl NCIntegrable for f64 {
    type NodeY = Vec<f64>;
    type NCPolynomial = Polynomial;

    fn compute_node_y<F>(f: F, node_x: &[f64]) -> Vec<f64>
    where
        F: Fn(f64) -> Self,
    {
        node_x.to_vec().fmap(|x| f(x))
    }

    fn compute_polynomial(node_x: &[f64], node_y: &Vec<f64>) -> Polynomial {
        lagrange_polynomial(node_x.to_vec(), node_y.to_vec())
    }

    fn integrate_polynomial(p: &Polynomial) -> Polynomial {
        p.integral()
    }

    fn evaluate_polynomial(p: &Polynomial, x: f64) -> Self {
        p.eval(x)
    }
}

/// Type that can be integrated using Gauss Legendre or Kronrod Quadrature
pub trait GLKIntegrable:
    std::ops::Add<Self, Output = Self> + std::ops::Mul<f64, Output = Self> + Sized
{
    /// The neutral element of addition.
    const ZERO: Self;
}

impl GLKIntegrable for f64 {
    const ZERO: Self = 0f64;
}

/// Type that can be integrated using Gauss Kronrod Quadrature
pub trait GKIntegrable: GLKIntegrable + std::ops::Sub<Self, Output = Self> + Sized + Clone {
    /// Returns `true` if this object is neither infinite nor NaN.
    fn is_finite(&self) -> bool;

    /// Returns the value of the norm used by Gauss Kronrod Quadrature to
    /// compute relative tolerances.
    fn gk_norm(&self) -> f64;

    /// Returns true if this object is less than the tolerance passed as
    /// the `tol` argument. By default, returns true if its [`gk_norm`](GKIntegrable::gk_norm)
    /// is less than `tol`.
    fn is_within_tol(&self, tol: f64) -> bool {
        self.gk_norm() < tol
    }
}

impl GKIntegrable for f64 {
    fn is_finite(&self) -> bool {
        f64::is_finite(*self)
    }

    fn gk_norm(&self) -> f64 {
        self.abs()
    }
}

/// Numerical Integration
///
/// # Description
/// `fn integrate(f, (a,b), method) -> Y`
///
/// `Y` must implement [`GKIntegrable`] and [`NCIntegrable`], like `f64` or
/// `C64` (`complex` feature only).
///
/// * `f`: Target function (`Fn(f64) -> Y`)
/// * `(a,b)` : Target interval
/// * `method` : Numerical integration method
///
/// # Method
///
/// * Gauss-Legendre Quadrature (up to order 30) : `GaussLegendre(usize)`
/// * Newton-Cotes Quadrature: `NewtonCotes(usize)`
/// * Gauss-Kronrod Quadrature
///     * `G7K15(tol, max_iter)`
///     * `G10K21`
///     * `G15K31`
///     * `G20K41`
///     * `G25K51`
///     * `G30K61`
/// * Gauss-Kronrod Quarature with Relative Error
///     * `G7K15R(rtol, max_iter)`
///     * `G10K21R`
///     * `G15K31R`
///     * `G20K41R`
///     * `G25K51R`
///     * `G30K61R`
pub fn integrate<F, Y>(f: F, (a, b): (f64, f64), method: Integral) -> Y
where
    F: Fn(f64) -> Y + Copy,
    Y: GKIntegrable + NCIntegrable,
{
    match method {
        Integral::GaussLegendre(n) => gauss_legendre_quadrature(f, n, (a, b)),
        Integral::NewtonCotes(n) => newton_cotes_quadrature(f, n, (a, b)),
        method => gauss_kronrod_quadrature_optimized(f, (a, b), method),
    }
}

/// Newton Cotes Quadrature
pub fn newton_cotes_quadrature<F, Y>(f: F, n: usize, (a, b): (f64, f64)) -> Y
where
    F: Fn(f64) -> Y,
    Y: NCIntegrable,
{
    let h = (b - a) / (n as f64);
    let node_x = seq(a, b, h);
    let node_y = Y::compute_node_y(f, &node_x);
    let p = Y::compute_polynomial(&node_x, &node_y);
    let q = Y::integrate_polynomial(&p);

    Y::evaluate_polynomial(&q, b) - Y::evaluate_polynomial(&q, a)
}

/// Gauss Legendre Quadrature
///
/// # Type
/// * `f, n, (a,b) -> Y`
///     * `Y` implements [`GLKIntegrable`], like `f64` or `C64` (`complex` feature
///       only)
///     * `f`: Numerical function (`Fn(f64) -> Y`)
///     * `n`: Order of Legendre polynomial (up to 16)
///     * `(a,b)`: Interval of integration
///
/// # Reference
/// * A. N. Lowan et al. (1942)
pub fn gauss_legendre_quadrature<F, Y>(f: F, n: usize, (a, b): (f64, f64)) -> Y
where
    F: Fn(f64) -> Y,
    Y: GLKIntegrable,
{
    unit_gauss_legendre_quadrature(|x| f(x * (b - a) / 2f64 + (a + b) / 2f64), n) * ((b - a) / 2f64)
}

/// Gauss Kronrod Quadrature
///
/// # Type
/// * `f, (a,b), method -> Y`
///     * `Y` implements [`GKIntegrable`], like `f64` or `C64` (`complex` feature
///       only)
///     * `f`: Numerical function (`Fn(f64) -> Y + Copy`)
///     * `(a,b)`: Interval of integration
///     * `method`: Integration method
///         * `G7K15(tol)`
///
/// # Reference
/// * arXiv: [1003.4629](https://arxiv.org/abs/1003.4629)
#[allow(non_snake_case)]
pub fn gauss_kronrod_quadrature<F, T, S, Y>(f: F, (a, b): (T, S), method: Integral) -> Y
where
    F: Fn(f64) -> Y + Copy,
    T: Into<f64>,
    S: Into<f64>,
    Y: GKIntegrable,
{
    let (g, k) = method.get_gauss_kronrod_order();
    let tol = method.get_tol();
    let max_iter = method.get_max_iter();
    let mut I = Y::ZERO;
    let mut S: Vec<(f64, f64, f64, u32)> = vec![];
    S.push((a.into(), b.into(), tol, max_iter));

    loop {
        match S.pop() {
            Some((a, b, tol, max_iter)) => {
                let G = gauss_legendre_quadrature(f, g as usize, (a, b));
                let K = kronrod_quadrature(f, k as usize, (a, b));
                let c = (a + b) / 2f64;
                let tol_curr = if method.is_relative() {
                    tol * G.gk_norm()
                } else {
                    tol
                };
                if (G.clone() - K).is_within_tol(tol_curr) || a == b || max_iter == 0 {
                    if !G.is_finite() {
                        return G;
                    }
                    I = I + G;
                } else {
                    S.push((a, c, tol / 2f64, max_iter - 1));
                    S.push((c, b, tol / 2f64, max_iter - 1));
                }
            }
            None => break,
        }
    }
    I
}

pub fn kronrod_quadrature<F, Y>(f: F, n: usize, (a, b): (f64, f64)) -> Y
where
    F: Fn(f64) -> Y,
    Y: GLKIntegrable,
{
    unit_kronrod_quadrature(|x| f(x * (b - a) / 2f64 + (a + b) / 2f64), n) * ((b - a) / 2f64)
}

#[allow(non_snake_case)]
pub fn gauss_kronrod_quadrature_optimized<F, T, S, Y>(f: F, (a, b): (T, S), method: Integral) -> Y
where
    F: Fn(f64) -> Y,
    T: Into<f64>,
    S: Into<f64>,
    Y: GKIntegrable, // Requires Clone, Add, Sub, Mul<f64>, ZERO, gk_norm, is_finite
{
    // Extract parameters from the method enum
    let (g, k) = method.get_gauss_kronrod_order();
    let tol = method.get_tol(); // Base tolerance
    let max_iter = method.get_max_iter(); // Max recursion depth
    let is_relative = method.is_relative(); // Use relative tolerance?

    let mut I = Y::ZERO; // Accumulator for the integral result
    let mut S: Vec<(f64, f64, f64, u32)> = Vec::with_capacity(max_iter as usize + 10); // Stack for subintervals: (a, b, tolerance, depth)

    let a_f64 = a.into();
    let b_f64 = b.into();

    // Handle zero-width interval initially
    if a_f64 == b_f64 {
        return Y::ZERO;
    }

    // Push the initial interval onto the stack
    S.push((a_f64, b_f64, tol, max_iter));

    // Adaptive quadrature loop
    while let Some((a_curr, b_curr, tol_curr, iter_left)) = S.pop() {
        // Call the helper function to get G, K, and error estimate for the current interval
        // This performs all necessary function evaluations efficiently
        let (G, K, error_estimate_raw) =
            compute_gauss_kronrod_sums_stored(&f, (a_curr, b_curr), g as usize, k as usize);

        // Calculate the norm of the error estimate |G - K|
        let error_norm = error_estimate_raw.gk_norm();

        // Determine the tolerance threshold for this interval
        let current_tolerance_threshold = if is_relative {
            // Scale relative tolerance by the norm of the higher-order estimate (K)
            tol_curr * K.clone().gk_norm()
        } else {
            // Use absolute tolerance directly
            tol_curr
        };

        // Check termination conditions for this interval:
        // 1. Error is within tolerance
        // 2. Interval width is zero (or numerically close)
        // 3. Maximum iteration depth reached
        if error_norm <= current_tolerance_threshold || a_curr == b_curr || iter_left == 0 {
            // Interval meets criteria, add its contribution to the total integral I
            // Prefer the higher-order Kronrod estimate (K) if it's finite.
            if !K.is_finite() {
                // If K is not finite, check G
                if G.is_finite() {
                    // If K is not finite but G is, add G
                    I = I + G.clone();
                }
            } else {
                // K is finite, add it to the total integral
                I = I + K.clone();
            }
        } else {
            // Interval does not meet criteria, subdivide it
            let c = (a_curr + b_curr) * 0.5; // Midpoint

            // Avoid infinite loops if midpoint calculation doesn't change due to precision
            if c <= a_curr || c >= b_curr {
                // Cannot subdivide further, precision limit reached.
                // Add the best estimate (K if finite, else G if finite) for this interval and stop subdividing.
                if !K.is_finite() {
                    if G.is_finite() {
                        I = I + G.clone();
                    }
                } else {
                    I = I + K.clone();
                }
                continue; // Skip pushing to stack
            }

            // Tolerance propagation for subintervals: divide by sqrt(2)
            // This is often preferred over dividing by 2 for better error distribution.
            let next_tol = tol_curr / 1.4142135623730951; // tol / sqrt(2)
            let next_iter = iter_left - 1; // Decrement remaining iterations

            // Push the two subintervals onto the stack for further processing
            S.push((a_curr, c, next_tol, next_iter)); // Left subinterval [a, c]
            S.push((c, b_curr, next_tol, next_iter)); // Right subinterval [c, b]
        }
    }

    I
}

/// Helper function to compute G and K integrals reusing function evaluations.
/// Optimized to store all function evaluations.
/// Returns (Gauss Integral G, Kronrod Integral K, Error Estimate G - K).
#[allow(non_snake_case)]
fn compute_gauss_kronrod_sums_stored<F, Y>(
    f: F,
    (a, b): (f64, f64),
    g: usize, // Gauss order
    k: usize, // Kronrod order
) -> (Y, Y, Y)
where
    F: Fn(f64) -> Y,
    Y: GKIntegrable,
{
    // 1. Get Kronrod nodes and weights
    // Assumes nodes ordered [-xn,...,0,...,+xn] and weights correspond.
    let (kronrod_weights, kronrod_nodes) = kronrod_table(k);
    // 2. Get Gauss weights
    // Assumes weights ordered [-wgn,...,wg0,...,+wgn] and symmetric w(+x) == w(-x).
    let (gauss_weights, _gauss_nodes_unused) = gauss_legendre_table(g);

    // 3. Calculate interval midpoint and half-length
    let xm = (a + b) * 0.5;
    let xh = (b - a) * 0.5;

    // Handle zero interval length
    if xh == 0.0 {
        let center_f = f(xm);
        return (center_f.clone(), center_f, Y::ZERO);
    }

    // 4. Evaluate f at all Kronrod nodes and store results (matching node order)
    let mut f_evals = Vec::with_capacity(k);
    for i in 0..k {
        let node = kronrod_nodes[i];
        f_evals.push(f(xm + xh * node));
    }

    // 5. Calculate Kronrod sum (K)
    let mut kronrod_sum = Y::ZERO;
    for i in 0..k {
        kronrod_sum = kronrod_sum + f_evals[i].clone() * kronrod_weights[i];
    }
    let K = kronrod_sum * xh;

    // 6. Calculate Gauss sum (G) using stored evaluations and Gauss weights
    // Assumes G(g) nodes are subset of K(k), g is odd, standard mapping: +/- xg_j maps to +/- xk_{2j}.
    let kronrod_center_idx = k / 2; // Index of 0.0 node in Kronrod/f_evals
    let gauss_center_idx = g / 2; // Index of 0.0 weight in Gauss weights
    let num_gauss_pairs = (g - 1) / 2;

    // Contribution from center node (0.0)
    let mut gauss_sum = f_evals[kronrod_center_idx].clone() * gauss_weights[gauss_center_idx];

    // Contribution from paired Gauss nodes (+/- xg_j)
    for j in 1..=num_gauss_pairs {
        if 2 * j > kronrod_center_idx {
            panic!(
                "Kronrod node index underflow. k={}, center_idx={}, j={}",
                k, kronrod_center_idx, j
            );
        }

        // Indices in f_evals corresponding to Kronrod nodes +/- xk_{2j}
        let kronrod_node_pos_idx = kronrod_center_idx + 2 * j;
        let kronrod_node_neg_idx = kronrod_center_idx - 2 * j;

        // Index in gauss_weights for the symmetric weight of the pair +/- xg_j
        let gauss_weight_pair_idx = gauss_center_idx + j;

        // Safety checks
        if kronrod_node_pos_idx >= k || kronrod_node_neg_idx >= k {
            panic!("Kronrod node index for Gauss sum out of bounds. Check mapping/indices. k={}, pos_idx={}, neg_idx={}", k, kronrod_node_pos_idx, kronrod_node_neg_idx);
        }
        if gauss_weight_pair_idx >= g {
            panic!("Gauss weight index out of bounds. Check gauss_weights structure. g={}, weight_idx={}", g, gauss_weight_pair_idx);
        }

        // Retrieve stored function evaluations
        let f_pos = f_evals[kronrod_node_pos_idx].clone();
        let f_neg = f_evals[kronrod_node_neg_idx].clone();

        // Retrieve the Gauss weight for the pair
        let weight_g_pair = gauss_weights[gauss_weight_pair_idx];

        // Add contribution: (f(+x) + f(-x)) * w_pair
        let f_sum_pair = f_pos + f_neg;
        gauss_sum = gauss_sum + f_sum_pair * weight_g_pair;
    }

    let G = gauss_sum * xh;

    // 7. Calculate Error Estimate (G - K)
    let error_estimate = G.clone() - K.clone();

    (G, K, error_estimate)
}

// =============================================================================
// Gauss Legendre Backends
// =============================================================================
fn unit_gauss_legendre_quadrature<F, Y>(f: F, n: usize) -> Y
where
    F: Fn(f64) -> Y,
    Y: GLKIntegrable,
{
    let (a, x) = gauss_legendre_table(n);
    let mut s = Y::ZERO;
    for i in 0..a.len() {
        s = s + f(x[i]) * a[i];
    }
    s
}

fn gauss_legendre_table(n: usize) -> (&'static [f64], &'static [f64]) {
    let root: &'static [f64] = match n {
        2 => &LEGENDRE_ROOT_2[..],
        3 => &LEGENDRE_ROOT_3[..],
        4 => &LEGENDRE_ROOT_4[..],
        5 => &LEGENDRE_ROOT_5[..],
        6 => &LEGENDRE_ROOT_6[..],
        7 => &LEGENDRE_ROOT_7[..],
        8 => &LEGENDRE_ROOT_8[..],
        9 => &LEGENDRE_ROOT_9[..],
        10 => &LEGENDRE_ROOT_10[..],
        11 => &LEGENDRE_ROOT_11[..],
        12 => &LEGENDRE_ROOT_12[..],
        13 => &LEGENDRE_ROOT_13[..],
        14 => &LEGENDRE_ROOT_14[..],
        15 => &LEGENDRE_ROOT_15[..],
        16 => &LEGENDRE_ROOT_16[..],
        17 => &LEGENDRE_ROOT_17[..],
        18 => &LEGENDRE_ROOT_18[..],
        19 => &LEGENDRE_ROOT_19[..],
        20 => &LEGENDRE_ROOT_20[..],
        21 => &LEGENDRE_ROOT_21[..],
        22 => &LEGENDRE_ROOT_22[..],
        23 => &LEGENDRE_ROOT_23[..],
        24 => &LEGENDRE_ROOT_24[..],
        25 => &LEGENDRE_ROOT_25[..],
        26 => &LEGENDRE_ROOT_26[..],
        27 => &LEGENDRE_ROOT_27[..],
        28 => &LEGENDRE_ROOT_28[..],
        29 => &LEGENDRE_ROOT_29[..],
        30 => &LEGENDRE_ROOT_30[..],
        _ => panic!("Legendre quadrature is limited up to n = 16"),
    };

    let weight: &'static [f64] = match n {
        2 => &LEGENDRE_WEIGHT_2[..],
        3 => &LEGENDRE_WEIGHT_3[..],
        4 => &LEGENDRE_WEIGHT_4[..],
        5 => &LEGENDRE_WEIGHT_5[..],
        6 => &LEGENDRE_WEIGHT_6[..],
        7 => &LEGENDRE_WEIGHT_7[..],
        8 => &LEGENDRE_WEIGHT_8[..],
        9 => &LEGENDRE_WEIGHT_9[..],
        10 => &LEGENDRE_WEIGHT_10[..],
        11 => &LEGENDRE_WEIGHT_11[..],
        12 => &LEGENDRE_WEIGHT_12[..],
        13 => &LEGENDRE_WEIGHT_13[..],
        14 => &LEGENDRE_WEIGHT_14[..],
        15 => &LEGENDRE_WEIGHT_15[..],
        16 => &LEGENDRE_WEIGHT_16[..],
        17 => &LEGENDRE_WEIGHT_17[..],
        18 => &LEGENDRE_WEIGHT_18[..],
        19 => &LEGENDRE_WEIGHT_19[..],
        20 => &LEGENDRE_WEIGHT_20[..],
        21 => &LEGENDRE_WEIGHT_21[..],
        22 => &LEGENDRE_WEIGHT_22[..],
        23 => &LEGENDRE_WEIGHT_23[..],
        24 => &LEGENDRE_WEIGHT_24[..],
        25 => &LEGENDRE_WEIGHT_25[..],
        26 => &LEGENDRE_WEIGHT_26[..],
        27 => &LEGENDRE_WEIGHT_27[..],
        28 => &LEGENDRE_WEIGHT_28[..],
        29 => &LEGENDRE_WEIGHT_29[..],
        30 => &LEGENDRE_WEIGHT_30[..],
        _ => panic!("Legendre quadrature is limited up to n = 16"),
    };

    (weight, root)
}

// =============================================================================
// Gauss Kronrod Backends
// =============================================================================

fn unit_kronrod_quadrature<F, Y>(f: F, n: usize) -> Y
where
    F: Fn(f64) -> Y,
    Y: GLKIntegrable,
{
    let (a, x) = kronrod_table(n);
    let mut s = Y::ZERO;
    for i in 0..a.len() {
        s = s + f(x[i]) * a[i];
    }
    s
}

fn kronrod_table(n: usize) -> (&'static [f64], &'static [f64]) {
    let node: &'static [f64] = match n {
        15 => &KRONROD_NODES_15[..],
        21 => &KRONROD_NODES_21[..],
        31 => &KRONROD_NODES_31[..],
        41 => &KRONROD_NODES_41[..],
        51 => &KRONROD_NODES_51[..],
        61 => &KRONROD_NODES_61[..],
        _ => panic!("Not yet implemented"),
    };
    let weight: &'static [f64] = match n {
        15 => &KRONROD_WEIGHTS_15[..],
        21 => &KRONROD_WEIGHTS_21[..],
        31 => &KRONROD_WEIGHTS_31[..],
        41 => &KRONROD_WEIGHTS_41[..],
        51 => &KRONROD_WEIGHTS_51[..],
        61 => &KRONROD_WEIGHTS_61[..],
        _ => panic!("Not yet implemented"),
    };

    (weight, node)
}

// =============================================================================
// Table for Gauss-Legendre Quadrature (ref. A. N. Lowan et al. (1942))
// =============================================================================
const LEGENDRE_ROOT_2: [f64; 2] = [-0.577350269189626, 0.577350269189626];
const LEGENDRE_ROOT_3: [f64; 3] = [-0.774596669241483, 0.0, 0.774596669241483];
const LEGENDRE_ROOT_4: [f64; 4] = [
    -0.861136311594053,
    -0.339981043584856,
    0.339981043584856,
    0.861136311594053,
];
const LEGENDRE_ROOT_5: [f64; 5] = [
    -0.906179845938664,
    -0.538469310105683,
    0.0,
    0.538469310105683,
    0.906179845938664,
];
const LEGENDRE_ROOT_6: [f64; 6] = [
    -0.932469514203152,
    -0.661209386466265,
    -0.238619186083197,
    0.238619186083197,
    0.661209386466265,
    0.932469514203152,
];
const LEGENDRE_ROOT_7: [f64; 7] = [
    -0.949107912342759,
    -0.741531185599394,
    -0.405845151377397,
    0.0,
    0.405845151377397,
    0.741531185599394,
    0.949107912342759,
];

const LEGENDRE_ROOT_8: [f64; 8] = [
    -0.960289856497536,
    -0.796666477413627,
    -0.525532409916329,
    -0.18343464249565,
    0.18343464249565,
    0.525532409916329,
    0.796666477413627,
    0.960289856497536,
];
const LEGENDRE_ROOT_9: [f64; 9] = [
    -0.968160239507626,
    -0.836031107326636,
    -0.61337143270059,
    -0.324253423403809,
    0.0,
    0.324253423403809,
    0.61337143270059,
    0.836031107326636,
    0.968160239507626,
];
const LEGENDRE_ROOT_10: [f64; 10] = [
    -0.973906528517172,
    -0.865063366688985,
    -0.679409568299024,
    -0.433395394129247,
    -0.148874338981631,
    0.148874338981631,
    0.433395394129247,
    0.679409568299024,
    0.865063366688985,
    0.973906528517172,
];
const LEGENDRE_ROOT_11: [f64; 11] = [
    -0.978228658146057,
    -0.887062599768095,
    -0.730152005574049,
    -0.519096129110681,
    -0.269543155952345,
    0.0,
    0.269543155952345,
    0.519096129110681,
    0.730152005574049,
    0.887062599768095,
    0.978228658146057,
];
const LEGENDRE_ROOT_12: [f64; 12] = [
    -0.981560634246719,
    -0.904117256370475,
    -0.769902674194305,
    -0.587317954286617,
    -0.36783149891818,
    -0.125333408511469,
    0.125333408511469,
    0.36783149891818,
    0.587317954286617,
    0.769902674194305,
    0.904117256370475,
    0.981560634246719,
];
const LEGENDRE_ROOT_13: [f64; 13] = [
    -0.984183054718588,
    -0.917598399222978,
    -0.80157809073331,
    -0.64234933944034,
    -0.448492751036447,
    -0.230458315955135,
    0.0,
    0.230458315955135,
    0.448492751036447,
    0.64234933944034,
    0.80157809073331,
    0.917598399222978,
    0.984183054718588,
];
const LEGENDRE_ROOT_14: [f64; 14] = [
    -0.986283808696812,
    -0.928434883663574,
    -0.827201315069765,
    -0.687292904811685,
    -0.515248636358154,
    -0.31911236892789,
    -0.108054948707344,
    0.108054948707344,
    0.31911236892789,
    0.515248636358154,
    0.687292904811685,
    0.827201315069765,
    0.928434883663574,
    0.986283808696812,
];
const LEGENDRE_ROOT_15: [f64; 15] = [
    -0.987992518020485,
    -0.937273392400706,
    -0.848206583410427,
    -0.72441773136017,
    -0.570972172608539,
    -0.394151347077563,
    -0.201194093997435,
    0.0,
    0.201194093997435,
    0.394151347077563,
    0.570972172608539,
    0.72441773136017,
    0.848206583410427,
    0.937273392400706,
    0.987992518020485,
];
const LEGENDRE_ROOT_16: [f64; 16] = [
    -0.98940093499165,
    -0.944575023073233,
    -0.865631202387832,
    -0.755404408355003,
    -0.617876244402644,
    -0.458016777657227,
    -0.281603550779259,
    -0.095012509837637,
    0.095012509837637,
    0.281603550779259,
    0.458016777657227,
    0.617876244402644,
    0.755404408355003,
    0.865631202387832,
    0.944575023073233,
    0.98940093499165,
];
const LEGENDRE_ROOT_17: [f64; 17] = [
    -0.9905754753144174,
    -0.9506755217687678,
    -0.8802391537269859,
    -0.7815140038968014,
    -0.6576711592166907,
    -0.5126905370864769,
    -0.3512317634538763,
    -0.17848418149584785,
    0.0,
    0.17848418149584785,
    0.3512317634538763,
    0.5126905370864769,
    0.6576711592166907,
    0.7815140038968014,
    0.8802391537269859,
    0.9506755217687678,
    0.9905754753144174,
];
const LEGENDRE_ROOT_18: [f64; 18] = [
    -0.9915651684209309,
    -0.9558239495713977,
    -0.8926024664975557,
    -0.8037049589725231,
    -0.6916870430603532,
    -0.5597708310739475,
    -0.41175116146284263,
    -0.2518862256915055,
    -0.0847750130417353,
    0.0847750130417353,
    0.2518862256915055,
    0.41175116146284263,
    0.5597708310739475,
    0.6916870430603532,
    0.8037049589725231,
    0.8926024664975557,
    0.9558239495713977,
    0.9915651684209309,
];
const LEGENDRE_ROOT_19: [f64; 19] = [
    -0.9924068438435844,
    -0.96020815213483,
    -0.9031559036148179,
    -0.8227146565371428,
    -0.7209661773352294,
    -0.600545304661681,
    -0.46457074137596094,
    -0.31656409996362983,
    -0.16035864564022537,
    0.0,
    0.16035864564022537,
    0.31656409996362983,
    0.46457074137596094,
    0.600545304661681,
    0.7209661773352294,
    0.8227146565371428,
    0.9031559036148179,
    0.96020815213483,
    0.9924068438435844,
];
const LEGENDRE_ROOT_20: [f64; 20] = [
    -0.9931285991850949,
    -0.9639719272779138,
    -0.912234428251326,
    -0.8391169718222188,
    -0.7463319064601508,
    -0.636053680726515,
    -0.5108670019508271,
    -0.37370608871541955,
    -0.22778585114164507,
    -0.07652652113349734,
    0.07652652113349734,
    0.22778585114164507,
    0.37370608871541955,
    0.5108670019508271,
    0.636053680726515,
    0.7463319064601508,
    0.8391169718222188,
    0.912234428251326,
    0.9639719272779138,
    0.9931285991850949,
];
const LEGENDRE_ROOT_21: [f64; 21] = [
    -0.9937521706203895,
    -0.9672268385663063,
    -0.9200993341504008,
    -0.8533633645833173,
    -0.7684399634756779,
    -0.6671388041974123,
    -0.5516188358872198,
    -0.4243421202074388,
    -0.2880213168024011,
    -0.1455618541608951,
    0.0,
    0.1455618541608951,
    0.2880213168024011,
    0.4243421202074388,
    0.5516188358872198,
    0.6671388041974123,
    0.7684399634756779,
    0.8533633645833173,
    0.9200993341504008,
    0.9672268385663063,
    0.9937521706203895,
];
const LEGENDRE_ROOT_22: [f64; 22] = [
    -0.9942945854823992,
    -0.9700604978354287,
    -0.926956772187174,
    -0.8658125777203002,
    -0.7878168059792081,
    -0.6944872631866827,
    -0.5876404035069116,
    -0.469355837986757,
    -0.34193582089208424,
    -0.20786042668822127,
    -0.06973927331972223,
    0.06973927331972223,
    0.20786042668822127,
    0.34193582089208424,
    0.469355837986757,
    0.5876404035069116,
    0.6944872631866827,
    0.7878168059792081,
    0.8658125777203002,
    0.926956772187174,
    0.9700604978354287,
    0.9942945854823992,
];
const LEGENDRE_ROOT_23: [f64; 23] = [
    -0.9947693349975522,
    -0.9725424712181152,
    -0.9329710868260161,
    -0.8767523582704416,
    -0.8048884016188399,
    -0.7186613631319502,
    -0.6196098757636461,
    -0.5095014778460075,
    -0.3903010380302908,
    -0.26413568097034495,
    -0.1332568242984661,
    0.0,
    0.1332568242984661,
    0.26413568097034495,
    0.3903010380302908,
    0.5095014778460075,
    0.6196098757636461,
    0.7186613631319502,
    0.8048884016188399,
    0.8767523582704416,
    0.9329710868260161,
    0.9725424712181152,
    0.9947693349975522,
];
const LEGENDRE_ROOT_24: [f64; 24] = [
    -0.9951872199970213,
    -0.9747285559713095,
    -0.9382745520027328,
    -0.8864155270044011,
    -0.820001985973903,
    -0.7401241915785544,
    -0.6480936519369755,
    -0.5454214713888396,
    -0.4337935076260451,
    -0.3150426796961634,
    -0.1911188674736163,
    -0.06405689286260563,
    0.06405689286260563,
    0.1911188674736163,
    0.3150426796961634,
    0.4337935076260451,
    0.5454214713888396,
    0.6480936519369755,
    0.7401241915785544,
    0.820001985973903,
    0.8864155270044011,
    0.9382745520027328,
    0.9747285559713095,
    0.9951872199970213,
];
const LEGENDRE_ROOT_25: [f64; 25] = [
    -0.9955569697904981,
    -0.9766639214595175,
    -0.9429745712289743,
    -0.8949919978782753,
    -0.833442628760834,
    -0.7592592630373576,
    -0.6735663684734684,
    -0.577662930241223,
    -0.473002731445715,
    -0.36117230580938786,
    -0.24386688372098844,
    -0.1228646926107104,
    0.0,
    0.1228646926107104,
    0.24386688372098844,
    0.36117230580938786,
    0.473002731445715,
    0.577662930241223,
    0.6735663684734684,
    0.7592592630373576,
    0.833442628760834,
    0.8949919978782753,
    0.9429745712289743,
    0.9766639214595175,
    0.9955569697904981,
];
const LEGENDRE_ROOT_26: [f64; 26] = [
    -0.9958857011456169,
    -0.978385445956471,
    -0.9471590666617142,
    -0.9026378619843071,
    -0.845445942788498,
    -0.7763859488206789,
    -0.6964272604199573,
    -0.6066922930176181,
    -0.5084407148245057,
    -0.4030517551234863,
    -0.2920048394859569,
    -0.17685882035689018,
    -0.05923009342931321,
    0.05923009342931321,
    0.17685882035689018,
    0.2920048394859569,
    0.4030517551234863,
    0.5084407148245057,
    0.6066922930176181,
    0.6964272604199573,
    0.7763859488206789,
    0.845445942788498,
    0.9026378619843071,
    0.9471590666617142,
    0.978385445956471,
    0.9958857011456169,
];
const LEGENDRE_ROOT_27: [f64; 27] = [
    -0.9961792628889886,
    -0.9799234759615012,
    -0.9509005578147051,
    -0.9094823206774911,
    -0.8562079080182945,
    -0.7917716390705082,
    -0.7170134737394237,
    -0.6329079719464952,
    -0.5405515645794569,
    -0.44114825175002687,
    -0.3359939036385089,
    -0.22645936543953685,
    -0.11397258560952997,
    0.0,
    0.11397258560952997,
    0.22645936543953685,
    0.3359939036385089,
    0.44114825175002687,
    0.5405515645794569,
    0.6329079719464952,
    0.7170134737394237,
    0.7917716390705082,
    0.8562079080182945,
    0.9094823206774911,
    0.9509005578147051,
    0.9799234759615012,
    0.9961792628889886,
];
const LEGENDRE_ROOT_28: [f64; 28] = [
    -0.9964424975739544,
    -0.9813031653708727,
    -0.9542592806289382,
    -0.9156330263921321,
    -0.8658925225743951,
    -0.8056413709171791,
    -0.7356108780136318,
    -0.656651094038865,
    -0.5697204718114017,
    -0.4758742249551183,
    -0.3762515160890787,
    -0.2720616276351781,
    -0.16456928213338076,
    -0.05507928988403427,
    0.05507928988403427,
    0.16456928213338076,
    0.2720616276351781,
    0.3762515160890787,
    0.4758742249551183,
    0.5697204718114017,
    0.656651094038865,
    0.7356108780136318,
    0.8056413709171791,
    0.8658925225743951,
    0.9156330263921321,
    0.9542592806289382,
    0.9813031653708727,
    0.9964424975739544,
];
const LEGENDRE_ROOT_29: [f64; 29] = [
    -0.9966794422605966,
    -0.9825455052614132,
    -0.9572855957780877,
    -0.9211802329530587,
    -0.8746378049201028,
    -0.8181854876152524,
    -0.7524628517344771,
    -0.6782145376026865,
    -0.5962817971382278,
    -0.5075929551242276,
    -0.41315288817400864,
    -0.31403163786763993,
    -0.21135228616600107,
    -0.10627823013267923,
    0.0,
    0.10627823013267923,
    0.21135228616600107,
    0.31403163786763993,
    0.41315288817400864,
    0.5075929551242276,
    0.5962817971382278,
    0.6782145376026865,
    0.7524628517344771,
    0.8181854876152524,
    0.8746378049201028,
    0.9211802329530587,
    0.9572855957780877,
    0.9825455052614132,
    0.9966794422605966,
];
const LEGENDRE_ROOT_30: [f64; 30] = [
    -0.9968934840746495,
    -0.9836681232797472,
    -0.9600218649683075,
    -0.9262000474292743,
    -0.8825605357920527,
    -0.8295657623827684,
    -0.7677774321048262,
    -0.6978504947933158,
    -0.6205261829892429,
    -0.5366241481420199,
    -0.44703376953808915,
    -0.3527047255308781,
    -0.25463692616788985,
    -0.15386991360858354,
    -0.0514718425553177,
    0.0514718425553177,
    0.15386991360858354,
    0.25463692616788985,
    0.3527047255308781,
    0.44703376953808915,
    0.5366241481420199,
    0.6205261829892429,
    0.6978504947933158,
    0.7677774321048262,
    0.8295657623827684,
    0.8825605357920527,
    0.9262000474292743,
    0.9600218649683075,
    0.9836681232797472,
    0.9968934840746495,
];

const LEGENDRE_WEIGHT_2: [f64; 2] = [1.0, 1.0];
const LEGENDRE_WEIGHT_3: [f64; 3] = [0.555555555555556, 0.888888888888889, 0.555555555555556];
const LEGENDRE_WEIGHT_4: [f64; 4] = [
    0.347854845137454,
    0.652145154862546,
    0.652145154862546,
    0.347854845137454,
];
const LEGENDRE_WEIGHT_5: [f64; 5] = [
    0.236926885056189,
    0.478628670499366,
    0.568888888888889,
    0.478628670499366,
    0.236926885056189,
];
const LEGENDRE_WEIGHT_6: [f64; 6] = [
    0.17132449237917,
    0.360761573048139,
    0.467913934572691,
    0.467913934572691,
    0.360761573048139,
    0.17132449237917,
];
const LEGENDRE_WEIGHT_7: [f64; 7] = [
    0.12948496616887,
    0.279705391489277,
    0.381830050505119,
    0.417959183673469,
    0.381830050505119,
    0.279705391489277,
    0.12948496616887,
];
const LEGENDRE_WEIGHT_8: [f64; 8] = [
    0.101228536290376,
    0.222381034453374,
    0.313706645877887,
    0.362683783378362,
    0.362683783378362,
    0.313706645877887,
    0.222381034453374,
    0.101228536290376,
];
const LEGENDRE_WEIGHT_9: [f64; 9] = [
    0.081274388361574,
    0.180648160694857,
    0.260610696402935,
    0.312347077040003,
    0.33023935500126,
    0.312347077040003,
    0.260610696402935,
    0.180648160694857,
    0.081274388361574,
];
const LEGENDRE_WEIGHT_10: [f64; 10] = [
    0.066671344308688,
    0.149451349150581,
    0.219086362515982,
    0.269266719309996,
    0.295524224714753,
    0.295524224714753,
    0.269266719309996,
    0.219086362515982,
    0.149451349150581,
    0.066671344308688,
];
const LEGENDRE_WEIGHT_11: [f64; 11] = [
    0.055668567116174,
    0.125580369464905,
    0.186290210927734,
    0.23319376459199,
    0.262804544510247,
    0.272925086777901,
    0.262804544510247,
    0.23319376459199,
    0.186290210927734,
    0.125580369464905,
    0.055668567116174,
];
const LEGENDRE_WEIGHT_12: [f64; 12] = [
    0.047175336386512,
    0.106939325995318,
    0.160078328543346,
    0.203167426723066,
    0.233492536538355,
    0.249147045813403,
    0.249147045813403,
    0.233492536538355,
    0.203167426723066,
    0.160078328543346,
    0.106939325995318,
    0.047175336386512,
];
const LEGENDRE_WEIGHT_13: [f64; 13] = [
    0.040484004765316,
    0.092121499837728,
    0.138873510219787,
    0.178145980761946,
    0.207816047536889,
    0.226283180262897,
    0.232551553230874,
    0.226283180262897,
    0.207816047536889,
    0.178145980761946,
    0.138873510219787,
    0.092121499837728,
    0.040484004765316,
];
const LEGENDRE_WEIGHT_14: [f64; 14] = [
    0.035119460331752,
    0.08015808715976,
    0.121518570687903,
    0.157203167158194,
    0.185538397477938,
    0.20519846372129,
    0.215263853463158,
    0.215263853463158,
    0.20519846372129,
    0.185538397477938,
    0.157203167158194,
    0.121518570687903,
    0.08015808715976,
    0.035119460331752,
];
const LEGENDRE_WEIGHT_15: [f64; 15] = [
    0.030753241996117,
    0.070366047488108,
    0.107159220467172,
    0.139570677926154,
    0.166269205816994,
    0.186161000015562,
    0.198431485327111,
    0.202578241925561,
    0.198431485327111,
    0.186161000015562,
    0.166269205816994,
    0.139570677926154,
    0.107159220467172,
    0.070366047488108,
    0.030753241996117,
];
const LEGENDRE_WEIGHT_16: [f64; 16] = [
    0.027152459411754,
    0.062253523938648,
    0.095158511682493,
    0.124628971255534,
    0.149595988816577,
    0.169156519395003,
    0.182603415044924,
    0.189450610455069,
    0.189450610455069,
    0.182603415044924,
    0.169156519395003,
    0.149595988816577,
    0.124628971255534,
    0.095158511682493,
    0.062253523938648,
    0.027152459411754,
];
const LEGENDRE_WEIGHT_17: [f64; 17] = [
    0.02414830286854793,
    0.0554595293739872,
    0.08503614831717918,
    0.11188384719340397,
    0.13513636846852548,
    0.15404576107681028,
    0.16800410215645004,
    0.17656270536699264,
    0.17944647035620653,
    0.17656270536699264,
    0.16800410215645004,
    0.15404576107681028,
    0.13513636846852548,
    0.11188384719340397,
    0.08503614831717918,
    0.0554595293739872,
    0.02414830286854793,
];
const LEGENDRE_WEIGHT_18: [f64; 18] = [
    0.02161601352648331,
    0.0497145488949698,
    0.07642573025488907,
    0.10094204410628717,
    0.12255520671147846,
    0.14064291467065065,
    0.15468467512626524,
    0.16427648374583273,
    0.1691423829631436,
    0.1691423829631436,
    0.16427648374583273,
    0.15468467512626524,
    0.14064291467065065,
    0.12255520671147846,
    0.10094204410628717,
    0.07642573025488907,
    0.0497145488949698,
    0.02161601352648331,
];
const LEGENDRE_WEIGHT_19: [f64; 19] = [
    0.019461788229726478,
    0.0448142267656996,
    0.06904454273764123,
    0.09149002162245,
    0.11156664554733399,
    0.12875396253933621,
    0.1426067021736066,
    0.15276604206585967,
    0.15896884339395434,
    0.1610544498487837,
    0.15896884339395434,
    0.15276604206585967,
    0.1426067021736066,
    0.12875396253933621,
    0.11156664554733399,
    0.09149002162245,
    0.06904454273764123,
    0.0448142267656996,
    0.019461788229726478,
];
const LEGENDRE_WEIGHT_20: [f64; 20] = [
    0.017614007139152118,
    0.04060142980038694,
    0.06267204833410907,
    0.08327674157670475,
    0.10193011981724044,
    0.11819453196151841,
    0.13168863844917664,
    0.14209610931838204,
    0.14917298647260374,
    0.15275338713072584,
    0.15275338713072584,
    0.14917298647260374,
    0.14209610931838204,
    0.13168863844917664,
    0.11819453196151841,
    0.10193011981724044,
    0.08327674157670475,
    0.06267204833410907,
    0.04060142980038694,
    0.017614007139152118,
];
const LEGENDRE_WEIGHT_21: [f64; 21] = [
    0.016017228257774335,
    0.036953789770852494,
    0.05713442542685721,
    0.0761001136283793,
    0.09344442345603386,
    0.10879729916714838,
    0.12183141605372853,
    0.13226893863333747,
    0.13988739479107315,
    0.14452440398997005,
    0.14608113364969041,
    0.14452440398997005,
    0.13988739479107315,
    0.13226893863333747,
    0.12183141605372853,
    0.10879729916714838,
    0.09344442345603386,
    0.0761001136283793,
    0.05713442542685721,
    0.036953789770852494,
    0.016017228257774335,
];
const LEGENDRE_WEIGHT_22: [f64; 22] = [
    0.014627995298272202,
    0.03377490158481416,
    0.052293335152683286,
    0.06979646842452049,
    0.08594160621706773,
    0.10041414444288096,
    0.11293229608053922,
    0.12325237681051242,
    0.13117350478706238,
    0.13654149834601517,
    0.13925187285563198,
    0.13925187285563198,
    0.13654149834601517,
    0.13117350478706238,
    0.12325237681051242,
    0.11293229608053922,
    0.10041414444288096,
    0.08594160621706773,
    0.06979646842452049,
    0.052293335152683286,
    0.03377490158481416,
    0.014627995298272202,
];
const LEGENDRE_WEIGHT_23: [f64; 23] = [
    0.013411859487141771,
    0.030988005856979445,
    0.04803767173108467,
    0.06423242140852585,
    0.07928141177671895,
    0.09291576606003515,
    0.10489209146454141,
    0.11499664022241136,
    0.12304908430672953,
    0.12890572218808216,
    0.1324620394046966,
    0.13365457218610619,
    0.1324620394046966,
    0.12890572218808216,
    0.12304908430672953,
    0.11499664022241136,
    0.10489209146454141,
    0.09291576606003515,
    0.07928141177671895,
    0.06423242140852585,
    0.04803767173108467,
    0.030988005856979445,
    0.013411859487141771,
];
const LEGENDRE_WEIGHT_24: [f64; 24] = [
    0.0123412297999872,
    0.028531388628933663,
    0.04427743881741981,
    0.05929858491543678,
    0.07334648141108031,
    0.08619016153195327,
    0.09761865210411388,
    0.10744427011596563,
    0.1155056680537256,
    0.12167047292780339,
    0.1258374563468283,
    0.12793819534675216,
    0.12793819534675216,
    0.1258374563468283,
    0.12167047292780339,
    0.1155056680537256,
    0.10744427011596563,
    0.09761865210411388,
    0.08619016153195327,
    0.07334648141108031,
    0.05929858491543678,
    0.04427743881741981,
    0.028531388628933663,
    0.0123412297999872,
];
const LEGENDRE_WEIGHT_25: [f64; 25] = [
    0.011393798501026288,
    0.026354986615032137,
    0.040939156701306316,
    0.054904695975835194,
    0.06803833381235691,
    0.08014070033500102,
    0.09102826198296365,
    0.10053594906705064,
    0.10851962447426365,
    0.11485825914571164,
    0.11945576353578477,
    0.12224244299031004,
    0.12317605372671545,
    0.12224244299031004,
    0.11945576353578477,
    0.11485825914571164,
    0.10851962447426365,
    0.10053594906705064,
    0.09102826198296365,
    0.08014070033500102,
    0.06803833381235691,
    0.054904695975835194,
    0.040939156701306316,
    0.026354986615032137,
    0.011393798501026288,
];
const LEGENDRE_WEIGHT_26: [f64; 26] = [
    0.01055137261734301,
    0.02441785109263191,
    0.03796238329436276,
    0.05097582529714781,
    0.06327404632957484,
    0.07468414976565975,
    0.08504589431348523,
    0.09421380035591415,
    0.10205916109442542,
    0.10847184052857659,
    0.11336181654631967,
    0.11666044348529658,
    0.11832141527926228,
    0.11832141527926228,
    0.11666044348529658,
    0.11336181654631967,
    0.10847184052857659,
    0.10205916109442542,
    0.09421380035591415,
    0.08504589431348523,
    0.07468414976565975,
    0.06327404632957484,
    0.05097582529714781,
    0.03796238329436276,
    0.02441785109263191,
    0.01055137261734301,
];
const LEGENDRE_WEIGHT_27: [f64; 27] = [
    0.00979899605129436,
    0.02268623159618062,
    0.03529705375741971,
    0.04744941252061506,
    0.0589835368598336,
    0.0697488237662456,
    0.07960486777305777,
    0.08842315854375694,
    0.0960887273700285,
    0.1025016378177458,
    0.10757828578853319,
    0.11125248835684519,
    0.11347634610896515,
    0.114220867378957,
    0.11347634610896515,
    0.11125248835684519,
    0.10757828578853319,
    0.1025016378177458,
    0.0960887273700285,
    0.08842315854375694,
    0.07960486777305777,
    0.0697488237662456,
    0.0589835368598336,
    0.04744941252061506,
    0.03529705375741971,
    0.02268623159618062,
    0.00979899605129436,
];
const LEGENDRE_WEIGHT_28: [f64; 28] = [
    0.00912428259309452,
    0.02113211259277126,
    0.03290142778230438,
    0.04427293475900423,
    0.05510734567571675,
    0.0652729239669996,
    0.07464621423456878,
    0.08311341722890121,
    0.09057174439303284,
    0.09693065799792992,
    0.10211296757806076,
    0.10605576592284642,
    0.10871119225829413,
    0.1100470130164752,
    0.1100470130164752,
    0.10871119225829413,
    0.10605576592284642,
    0.10211296757806076,
    0.09693065799792992,
    0.09057174439303284,
    0.08311341722890121,
    0.07464621423456878,
    0.0652729239669996,
    0.05510734567571675,
    0.04427293475900423,
    0.03290142778230438,
    0.02113211259277126,
    0.00912428259309452,
];
const LEGENDRE_WEIGHT_29: [f64; 29] = [
    0.00851690387874641,
    0.019732085056122707,
    0.030740492202093624,
    0.041402062518682836,
    0.05159482690249792,
    0.061203090657079136,
    0.07011793325505128,
    0.07823832713576379,
    0.08547225736617253,
    0.09173775713925876,
    0.0969638340944086,
    0.10109127375991497,
    0.10407331007772938,
    0.10587615509732094,
    0.10647938171831424,
    0.10587615509732094,
    0.10407331007772938,
    0.10109127375991497,
    0.0969638340944086,
    0.09173775713925876,
    0.08547225736617253,
    0.07823832713576379,
    0.07011793325505128,
    0.061203090657079136,
    0.05159482690249792,
    0.041402062518682836,
    0.030740492202093624,
    0.019732085056122707,
    0.00851690387874641,
];
const LEGENDRE_WEIGHT_30: [f64; 30] = [
    0.007968192496166607,
    0.01846646831109096,
    0.02878470788332337,
    0.03879919256962705,
    0.04840267283059405,
    0.057493156217619065,
    0.0659742298821805,
    0.0737559747377052,
    0.08075589522942021,
    0.08689978720108298,
    0.09212252223778614,
    0.09636873717464425,
    0.09959342058679527,
    0.1017623897484055,
    0.10285265289355884,
    0.10285265289355884,
    0.1017623897484055,
    0.09959342058679527,
    0.09636873717464425,
    0.09212252223778614,
    0.08689978720108298,
    0.08075589522942021,
    0.0737559747377052,
    0.0659742298821805,
    0.057493156217619065,
    0.04840267283059405,
    0.03879919256962705,
    0.02878470788332337,
    0.01846646831109096,
    0.007968192496166607,
];

// =============================================================================
// Table for Kronrod Quadrature
// =============================================================================
const KRONROD_NODES_15: [f64; 15] = [
    -0.9914553711208126,
    -0.9491079123427585,
    -0.8648644233597691,
    -0.7415311855993945,
    -0.5860872354676911,
    -0.4058451513773972,
    -0.20778495500789848,
    0.0,
    0.20778495500789848,
    0.4058451513773972,
    0.5860872354676911,
    0.7415311855993945,
    0.8648644233597691,
    0.9491079123427585,
    0.9914553711208126,
];

const KRONROD_WEIGHTS_15: [f64; 15] = [
    0.022935322010529224,
    0.06309209262997854,
    0.10479001032225017,
    0.14065325971552592,
    0.1690047266392679,
    0.19035057806478542,
    0.20443294007529889,
    0.20948214108472782,
    0.20443294007529889,
    0.19035057806478542,
    0.1690047266392679,
    0.14065325971552592,
    0.10479001032225017,
    0.06309209262997854,
    0.022935322010529224,
];

const KRONROD_NODES_21: [f64; 21] = [
    -0.9956571630258081,
    -0.9739065285171717,
    -0.9301574913557082,
    -0.8650633666889845,
    -0.7808177265864169,
    -0.6794095682990244,
    -0.5627571346686047,
    -0.4333953941292472,
    -0.2943928627014602,
    -0.14887433898163122,
    0.0,
    0.14887433898163122,
    0.2943928627014602,
    0.4333953941292472,
    0.5627571346686047,
    0.6794095682990244,
    0.7808177265864169,
    0.8650633666889845,
    0.9301574913557082,
    0.9739065285171717,
    0.9956571630258081,
];

const KRONROD_WEIGHTS_21: [f64; 21] = [
    0.011694638867371874,
    0.032558162307964725,
    0.054755896574351995,
    0.07503967481091996,
    0.09312545458369761,
    0.10938715880229764,
    0.12349197626206584,
    0.13470921731147334,
    0.14277593857706009,
    0.14773910490133849,
    0.1494455540029169,
    0.14773910490133849,
    0.14277593857706009,
    0.13470921731147334,
    0.12349197626206584,
    0.10938715880229764,
    0.09312545458369761,
    0.07503967481091996,
    0.054755896574351995,
    0.032558162307964725,
    0.011694638867371874,
];

const KRONROD_NODES_31: [f64; 31] = [
    -0.9980022986933971,
    -0.9879925180204854,
    -0.9677390756791391,
    -0.937273392400706,
    -0.8972645323440819,
    -0.8482065834104272,
    -0.790418501442466,
    -0.7244177313601701,
    -0.650996741297417,
    -0.5709721726085388,
    -0.4850818636402397,
    -0.3941513470775634,
    -0.29918000715316884,
    -0.20119409399743451,
    -0.1011420669187175,
    0.0,
    0.1011420669187175,
    0.20119409399743451,
    0.29918000715316884,
    0.3941513470775634,
    0.4850818636402397,
    0.5709721726085388,
    0.650996741297417,
    0.7244177313601701,
    0.790418501442466,
    0.8482065834104272,
    0.8972645323440819,
    0.937273392400706,
    0.9677390756791391,
    0.9879925180204854,
    0.9980022986933971,
];

const KRONROD_WEIGHTS_31: [f64; 31] = [
    0.005377479872923349,
    0.015007947329316124,
    0.02546084732671532,
    0.03534636079137585,
    0.04458975132476488,
    0.05348152469092809,
    0.06200956780067064,
    0.06985412131872826,
    0.07684968075772038,
    0.08308050282313302,
    0.08856444305621176,
    0.09312659817082532,
    0.09664272698362368,
    0.09917359872179196,
    0.10076984552387559,
    0.10133000701479154,
    0.10076984552387559,
    0.09917359872179196,
    0.09664272698362368,
    0.09312659817082532,
    0.08856444305621176,
    0.08308050282313302,
    0.07684968075772038,
    0.06985412131872826,
    0.06200956780067064,
    0.05348152469092809,
    0.04458975132476488,
    0.03534636079137585,
    0.02546084732671532,
    0.015007947329316124,
    0.005377479872923349,
];

const KRONROD_NODES_41: [f64; 41] = [
    -0.9988590315882777,
    -0.9931285991850949,
    -0.9815078774502503,
    -0.9639719272779138,
    -0.9408226338317548,
    -0.912234428251326,
    -0.878276811252282,
    -0.8391169718222188,
    -0.7950414288375512,
    -0.7463319064601508,
    -0.6932376563347514,
    -0.636053680726515,
    -0.5751404468197103,
    -0.5108670019508271,
    -0.4435931752387251,
    -0.37370608871541955,
    -0.301627868114913,
    -0.22778585114164507,
    -0.15260546524092267,
    -0.07652652113349734,
    0.0,
    0.07652652113349734,
    0.15260546524092267,
    0.22778585114164507,
    0.301627868114913,
    0.37370608871541955,
    0.4435931752387251,
    0.5108670019508271,
    0.5751404468197103,
    0.636053680726515,
    0.6932376563347514,
    0.7463319064601508,
    0.7950414288375512,
    0.8391169718222188,
    0.878276811252282,
    0.912234428251326,
    0.9408226338317548,
    0.9639719272779138,
    0.9815078774502503,
    0.9931285991850949,
    0.9988590315882777,
];

const KRONROD_WEIGHTS_41: [f64; 41] = [
    0.00307358371852053,
    0.008600269855642943,
    0.014626169256971253,
    0.020388373461266523,
    0.02588213360495116,
    0.0312873067770328,
    0.036600169758200796,
    0.04166887332797369,
    0.04643482186749767,
    0.05094457392372869,
    0.05519510534828599,
    0.05911140088063957,
    0.06265323755478117,
    0.06583459713361842,
    0.06864867292852161,
    0.07105442355344407,
    0.07303069033278667,
    0.07458287540049918,
    0.07570449768455667,
    0.07637786767208074,
    0.07660071191799966,
    0.07637786767208074,
    0.07570449768455667,
    0.07458287540049918,
    0.07303069033278667,
    0.07105442355344407,
    0.06864867292852161,
    0.06583459713361842,
    0.06265323755478117,
    0.05911140088063957,
    0.05519510534828599,
    0.05094457392372869,
    0.04643482186749767,
    0.04166887332797369,
    0.036600169758200796,
    0.0312873067770328,
    0.02588213360495116,
    0.020388373461266523,
    0.014626169256971253,
    0.008600269855642943,
    0.00307358371852053,
];

const KRONROD_NODES_51: [f64; 51] = [
    -0.9992621049926098,
    -0.9955569697904981,
    -0.9880357945340772,
    -0.9766639214595175,
    -0.9616149864258425,
    -0.9429745712289743,
    -0.9207471152817016,
    -0.8949919978782753,
    -0.8658470652932756,
    -0.833442628760834,
    -0.7978737979985001,
    -0.7592592630373576,
    -0.7177664068130843,
    -0.6735663684734684,
    -0.6268100990103174,
    -0.577662930241223,
    -0.5263252843347191,
    -0.473002731445715,
    -0.4178853821930377,
    -0.36117230580938786,
    -0.30308953893110785,
    -0.24386688372098844,
    -0.1837189394210489,
    -0.1228646926107104,
    -0.06154448300568508,
    0.0,
    0.06154448300568508,
    0.1228646926107104,
    0.1837189394210489,
    0.24386688372098844,
    0.30308953893110785,
    0.36117230580938786,
    0.4178853821930377,
    0.473002731445715,
    0.5263252843347191,
    0.577662930241223,
    0.6268100990103174,
    0.6735663684734684,
    0.7177664068130843,
    0.7592592630373576,
    0.7978737979985001,
    0.833442628760834,
    0.8658470652932756,
    0.8949919978782753,
    0.9207471152817016,
    0.9429745712289743,
    0.9616149864258425,
    0.9766639214595175,
    0.9880357945340772,
    0.9955569697904981,
    0.9992621049926098,
];

const KRONROD_WEIGHTS_51: [f64; 51] = [
    0.001987383892330316,
    0.005561932135356714,
    0.009473973386174152,
    0.013236229195571676,
    0.0168478177091283,
    0.02043537114588284,
    0.02400994560695322,
    0.02747531758785174,
    0.030792300167387487,
    0.03400213027432934,
    0.03711627148341554,
    0.04008382550403238,
    0.04287284502017005,
    0.04550291304992179,
    0.04798253713883671,
    0.05027767908071567,
    0.05236288580640747,
    0.05425112988854549,
    0.055950811220412316,
    0.057437116361567835,
    0.058689680022394206,
    0.05972034032417406,
    0.06053945537604586,
    0.06112850971705305,
    0.061471189871425316,
    0.061580818067832936,
    0.061471189871425316,
    0.06112850971705305,
    0.06053945537604586,
    0.05972034032417406,
    0.058689680022394206,
    0.057437116361567835,
    0.055950811220412316,
    0.05425112988854549,
    0.05236288580640747,
    0.05027767908071567,
    0.04798253713883671,
    0.04550291304992179,
    0.04287284502017005,
    0.04008382550403238,
    0.03711627148341554,
    0.03400213027432934,
    0.030792300167387487,
    0.02747531758785174,
    0.02400994560695322,
    0.02043537114588284,
    0.0168478177091283,
    0.013236229195571676,
    0.009473973386174152,
    0.005561932135356714,
    0.001987383892330316,
];

const KRONROD_NODES_61: [f64; 61] = [
    -0.9994844100504906,
    -0.9968934840746495,
    -0.9916309968704046,
    -0.9836681232797472,
    -0.9731163225011262,
    -0.9600218649683075,
    -0.94437444474856,
    -0.9262000474292743,
    -0.9055733076999078,
    -0.8825605357920527,
    -0.8572052335460612,
    -0.8295657623827684,
    -0.799727835821839,
    -0.7677774321048262,
    -0.7337900624532268,
    -0.6978504947933158,
    -0.6600610641266269,
    -0.6205261829892429,
    -0.5793452358263617,
    -0.5366241481420199,
    -0.49248046786177857,
    -0.44703376953808915,
    -0.4004012548303944,
    -0.3527047255308781,
    -0.30407320227362505,
    -0.25463692616788985,
    -0.20452511668230988,
    -0.15386991360858354,
    -0.10280693796673702,
    -0.0514718425553177,
    0.0,
    0.0514718425553177,
    0.10280693796673702,
    0.15386991360858354,
    0.20452511668230988,
    0.25463692616788985,
    0.30407320227362505,
    0.3527047255308781,
    0.4004012548303944,
    0.44703376953808915,
    0.49248046786177857,
    0.5366241481420199,
    0.5793452358263617,
    0.6205261829892429,
    0.6600610641266269,
    0.6978504947933158,
    0.7337900624532268,
    0.7677774321048262,
    0.799727835821839,
    0.8295657623827684,
    0.8572052335460612,
    0.8825605357920527,
    0.9055733076999078,
    0.9262000474292743,
    0.94437444474856,
    0.9600218649683075,
    0.9731163225011262,
    0.9836681232797472,
    0.9916309968704046,
    0.9968934840746495,
    0.9994844100504906,
];

const KRONROD_WEIGHTS_61: [f64; 61] = [
    0.001389013698677008,
    0.003890461127099884,
    0.006630703915931292,
    0.009273279659517762,
    0.011823015253496341,
    0.0143697295070458,
    0.016920889189053274,
    0.019414141193942382,
    0.02182803582160919,
    0.0241911620780806,
    0.0265099548823331,
    0.02875404876504129,
    0.030907257562387762,
    0.03298144705748373,
    0.034979338028060025,
    0.03688236465182123,
    0.038678945624727595,
    0.04037453895153596,
    0.041969810215164244,
    0.04345253970135607,
    0.04481480013316266,
    0.04605923827100699,
    0.04718554656929915,
    0.048185861757087126,
    0.04905543455502978,
    0.0497956834270742,
    0.05040592140278235,
    0.0508817958987496,
    0.051221547849258774,
    0.05142612853745902,
    0.05149472942945157,
    0.05142612853745902,
    0.051221547849258774,
    0.0508817958987496,
    0.05040592140278235,
    0.0497956834270742,
    0.04905543455502978,
    0.048185861757087126,
    0.04718554656929915,
    0.04605923827100699,
    0.04481480013316266,
    0.04345253970135607,
    0.041969810215164244,
    0.04037453895153596,
    0.038678945624727595,
    0.03688236465182123,
    0.034979338028060025,
    0.03298144705748373,
    0.030907257562387762,
    0.02875404876504129,
    0.0265099548823331,
    0.0241911620780806,
    0.02182803582160919,
    0.019414141193942382,
    0.016920889189053274,
    0.0143697295070458,
    0.011823015253496341,
    0.009273279659517762,
    0.006630703915931292,
    0.003890461127099884,
    0.001389013698677008,
];
