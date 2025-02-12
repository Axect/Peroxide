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
        method => gauss_kronrod_quadrature(f, (a, b), method),
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

fn gauss_legendre_table(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut result_root = vec![0f64; n];
    let mut result_weight = vec![0f64; n];
    let ref_root: &[f64] = match n {
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

    let ref_weight: &[f64] = match n {
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

    match n % 2 {
        0 => {
            for i in 0..ref_root.len() {
                result_root[i] = ref_root[i];
                result_weight[i] = ref_weight[i];
            }

            for i in ref_root.len()..n {
                result_root[i] = -ref_root[n - i - 1];
                result_weight[i] = ref_weight[n - i - 1];
            }
        }
        1 => {
            for i in 0..ref_root.len() {
                result_root[i] = ref_root[i];
                result_weight[i] = ref_weight[i];
            }

            for i in ref_root.len()..n {
                result_root[i] = -ref_root[n - i];
                result_weight[i] = ref_weight[n - i];
            }
        }
        _ => unreachable!(),
    }
    (result_weight, result_root)
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

fn kronrod_table(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut result_node = vec![0f64; n];
    let mut result_weight = vec![0f64; n];
    let ref_node: &[f64] = match n {
        15 => &KRONROD_NODES_15[..],
        21 => &KRONROD_NODES_21[..],
        31 => &KRONROD_NODES_31[..],
        41 => &KRONROD_NODES_41[..],
        51 => &KRONROD_NODES_51[..],
        61 => &KRONROD_NODES_61[..],
        _ => panic!("Not yet implemented"),
    };
    let ref_weight: &[f64] = match n {
        15 => &KRONROD_WEIGHTS_15[..],
        21 => &KRONROD_WEIGHTS_21[..],
        31 => &KRONROD_WEIGHTS_31[..],
        41 => &KRONROD_WEIGHTS_41[..],
        51 => &KRONROD_WEIGHTS_51[..],
        61 => &KRONROD_WEIGHTS_61[..],
        _ => panic!("Not yet implemented"),
    };

    match n % 2 {
        0 => {
            for i in 0..ref_node.len() {
                result_node[i] = ref_node[i];
                result_weight[i] = ref_weight[i];
            }

            for i in ref_node.len()..n {
                result_node[i] = -ref_node[n - i - 1];
                result_weight[i] = ref_weight[n - i - 1];
            }
        }
        1 => {
            for i in 0..ref_node.len() {
                result_node[i] = ref_node[i];
                result_weight[i] = ref_weight[i];
            }

            for i in ref_node.len()..n {
                result_node[i] = -ref_node[n - i];
                result_weight[i] = ref_weight[n - i];
            }
        }
        _ => unreachable!(),
    }
    (result_weight, result_node)
}

// =============================================================================
// Table for Gauss-Legendre Quadrature (ref. A. N. Lowan et al. (1942))
// =============================================================================
const LEGENDRE_ROOT_2: [f64; 1] = [0.577350269189626];
const LEGENDRE_ROOT_3: [f64; 2] = [0f64, 0.774596669241483];
const LEGENDRE_ROOT_4: [f64; 2] = [0.339981043584856, 0.861136311594053];
const LEGENDRE_ROOT_5: [f64; 3] = [0f64, 0.538469310105683, 0.906179845938664];
const LEGENDRE_ROOT_6: [f64; 3] = [0.238619186083197, 0.661209386466265, 0.932469514203152];
const LEGENDRE_ROOT_7: [f64; 4] = [
    0f64,
    0.405845151377397,
    0.741531185599394,
    0.949107912342759,
];
const LEGENDRE_ROOT_8: [f64; 4] = [
    0.183434642495650,
    0.525532409916329,
    0.796666477413627,
    0.960289856497536,
];
const LEGENDRE_ROOT_9: [f64; 5] = [
    0f64,
    0.324253423403809,
    0.613371432700590,
    0.836031107326636,
    0.968160239507626,
];
const LEGENDRE_ROOT_10: [f64; 5] = [
    0.148874338981631,
    0.433395394129247,
    0.679409568299024,
    0.865063366688985,
    0.973906528517172,
];
const LEGENDRE_ROOT_11: [f64; 6] = [
    0f64,
    0.269543155952345,
    0.519096129110681,
    0.730152005574049,
    0.887062599768095,
    0.978228658146057,
];
const LEGENDRE_ROOT_12: [f64; 6] = [
    0.125333408511469,
    0.367831498918180,
    0.587317954286617,
    0.769902674194305,
    0.904117256370475,
    0.981560634246719,
];
const LEGENDRE_ROOT_13: [f64; 7] = [
    0f64,
    0.230458315955135,
    0.448492751036447,
    0.642349339440340,
    0.801578090733310,
    0.917598399222978,
    0.984183054718588,
];
const LEGENDRE_ROOT_14: [f64; 7] = [
    0.108054948707344,
    0.319112368927890,
    0.515248636358154,
    0.687292904811685,
    0.827201315069765,
    0.928434883663574,
    0.986283808696812,
];
const LEGENDRE_ROOT_15: [f64; 8] = [
    0f64,
    0.201194093997435,
    0.394151347077563,
    0.570972172608539,
    0.724417731360170,
    0.848206583410427,
    0.937273392400706,
    0.987992518020485,
];
const LEGENDRE_ROOT_16: [f64; 8] = [
    0.095012509837637,
    0.281603550779259,
    0.458016777657227,
    0.617876244402644,
    0.755404408355003,
    0.865631202387832,
    0.944575023073233,
    0.989400934991650,
];
const LEGENDRE_ROOT_17: [f64; 9] = [
    0f64,
    0.178484181495847856,
    0.351231763453876315,
    0.512690537086476968,
    0.657671159216690766,
    0.781514003896801407,
    0.880239153726985902,
    0.950675521768767761,
    0.99057547531441734,
];
const LEGENDRE_ROOT_18: [f64; 9] = [
    0.084775013041735301,
    0.25188622569150551,
    0.411751161462842646,
    0.559770831073947535,
    0.691687043060353208,
    0.803704958972523116,
    0.892602466497555739,
    0.955823949571397755,
    0.991565168420930947,
];
const LEGENDRE_ROOT_19: [f64; 10] = [
    0f64,
    0.160358645640225376,
    0.316564099963629832,
    0.464570741375960946,
    0.600545304661681024,
    0.720966177335229379,
    0.822714656537142825,
    0.903155903614817902,
    0.960208152134830031,
    0.992406843843584403,
];
const LEGENDRE_ROOT_20: [f64; 10] = [
    0.076526521133497334,
    0.227785851141645078,
    0.373706088715419561,
    0.510867001950827098,
    0.636053680726515026,
    0.746331906460150793,
    0.839116971822218823,
    0.912234428251325906,
    0.963971927277913791,
    0.993128599185094925,
];
const LEGENDRE_ROOT_21: [f64; 11] = [
    0f64,
    0.145561854160895091,
    0.288021316802401097,
    0.424342120207438784,
    0.551618835887219807,
    0.667138804197412319,
    0.768439963475677909,
    0.853363364583317284,
    0.920099334150400829,
    0.967226838566306294,
    0.9937521706203895,
];
const LEGENDRE_ROOT_22: [f64; 11] = [
    0.069739273319722221,
    0.207860426688221285,
    0.341935820892084225,
    0.469355837986757026,
    0.58764040350691159,
    0.69448726318668278,
    0.787816805979208162,
    0.865812577720300137,
    0.926956772187174001,
    0.970060497835428727,
    0.994294585482399292,
];
const LEGENDRE_ROOT_23: [f64; 12] = [
    0f64,
    0.133256824298466111,
    0.26413568097034493,
    0.390301038030290831,
    0.50950147784600755,
    0.619609875763646156,
    0.718661363131950195,
    0.804888401618839892,
    0.876752358270441667,
    0.932971086826016102,
    0.972542471218115232,
    0.994769334997552124,
];
const LEGENDRE_ROOT_24: [f64; 12] = [
    0.064056892862605626,
    0.191118867473616309,
    0.315042679696163374,
    0.433793507626045138,
    0.545421471388839536,
    0.648093651936975569,
    0.740124191578554364,
    0.820001985973902922,
    0.886415527004401034,
    0.938274552002732759,
    0.974728555971309498,
    0.99518721999702136,
];
const LEGENDRE_ROOT_25: [f64; 13] = [
    0f64,
    0.122864692610710396,
    0.243866883720988432,
    0.361172305809387838,
    0.47300273144571496,
    0.577662930241222968,
    0.673566368473468365,
    0.759259263037357631,
    0.833442628760834001,
    0.894991997878275369,
    0.942974571228974339,
    0.976663921459517512,
    0.995556969790498098,
];
const LEGENDRE_ROOT_26: [f64; 13] = [
    0.059230093429313207,
    0.176858820356890184,
    0.292004839485956895,
    0.40305175512348631,
    0.508440714824505718,
    0.606692293017618063,
    0.696427260419957265,
    0.776385948820678856,
    0.845445942788498019,
    0.902637861984307074,
    0.94715906666171425,
    0.97838544595647099,
    0.99588570114561693,
];
const LEGENDRE_ROOT_27: [f64; 14] = [
    0f64,
    0.113972585609529967,
    0.226459365439536859,
    0.3359939036385089,
    0.441148251750026881,
    0.54055156457945689,
    0.632907971946495141,
    0.717013473739423699,
    0.791771639070508227,
    0.85620790801829449,
    0.909482320677491104,
    0.950900557814705007,
    0.979923475961501223,
    0.996179262888988567,
];
const LEGENDRE_ROOT_28: [f64; 14] = [
    0.0550792898840342704,
    0.164569282133380771,
    0.272061627635178078,
    0.37625151608907871,
    0.475874224955118261,
    0.569720471811401719,
    0.656651094038864961,
    0.735610878013631772,
    0.805641370917179171,
    0.865892522574395049,
    0.915633026392132074,
    0.954259280628938197,
    0.981303165370872754,
    0.99644249757395445,
];
const LEGENDRE_ROOT_29: [f64; 15] = [
    0f64,
    0.10627823013267923,
    0.211352286166001075,
    0.314031637867639935,
    0.413152888174008664,
    0.507592955124227642,
    0.59628179713822782,
    0.678214537602686515,
    0.752462851734477134,
    0.81818548761525245,
    0.87463780492010279,
    0.921180232953058785,
    0.957285595778087726,
    0.982545505261413175,
    0.996679442260596586,
];
const LEGENDRE_ROOT_30: [f64; 15] = [
    0.051471842555317696,
    0.153869913608583547,
    0.254636926167889846,
    0.352704725530878113,
    0.447033769538089177,
    0.536624148142019899,
    0.620526182989242861,
    0.6978504947933158,
    0.767777432104826195,
    0.829565762382768397,
    0.882560535792052682,
    0.926200047429274326,
    0.960021864968307512,
    0.98366812327974721,
    0.99689348407464954,
];
const LEGENDRE_WEIGHT_2: [f64; 1] = [1f64];
const LEGENDRE_WEIGHT_3: [f64; 2] = [0.888888888888889, 0.555555555555556];
const LEGENDRE_WEIGHT_4: [f64; 2] = [0.652145154862546, 0.347854845137454];
const LEGENDRE_WEIGHT_5: [f64; 3] = [0.568888888888889, 0.478628670499366, 0.236926885056189];
const LEGENDRE_WEIGHT_6: [f64; 3] = [0.467913934572691, 0.360761573048139, 0.171324492379170];
const LEGENDRE_WEIGHT_7: [f64; 4] = [
    0.417959183673469,
    0.381830050505119,
    0.279705391489277,
    0.129484966168870,
];
const LEGENDRE_WEIGHT_8: [f64; 4] = [
    0.362683783378362,
    0.313706645877887,
    0.222381034453374,
    0.101228536290376,
];
const LEGENDRE_WEIGHT_9: [f64; 5] = [
    0.330239355001260,
    0.312347077040003,
    0.260610696402935,
    0.180648160694857,
    0.081274388361574,
];
const LEGENDRE_WEIGHT_10: [f64; 5] = [
    0.295524224714753,
    0.269266719309996,
    0.219086362515982,
    0.149451349150581,
    0.066671344308688,
];
const LEGENDRE_WEIGHT_11: [f64; 6] = [
    0.272925086777901,
    0.262804544510247,
    0.233193764591990,
    0.186290210927734,
    0.125580369464905,
    0.055668567116174,
];
const LEGENDRE_WEIGHT_12: [f64; 6] = [
    0.249147045813403,
    0.233492536538355,
    0.203167426723066,
    0.160078328543346,
    0.106939325995318,
    0.047175336386512,
];
const LEGENDRE_WEIGHT_13: [f64; 7] = [
    0.232551553230874,
    0.226283180262897,
    0.207816047536889,
    0.178145980761946,
    0.138873510219787,
    0.092121499837728,
    0.040484004765316,
];
const LEGENDRE_WEIGHT_14: [f64; 7] = [
    0.215263853463158,
    0.205198463721290,
    0.185538397477938,
    0.157203167158194,
    0.121518570687903,
    0.080158087159760,
    0.035119460331752,
];
const LEGENDRE_WEIGHT_15: [f64; 8] = [
    0.202578241925561,
    0.198431485327111,
    0.186161000015562,
    0.166269205816994,
    0.139570677926154,
    0.107159220467172,
    0.070366047488108,
    0.030753241996117,
];
const LEGENDRE_WEIGHT_16: [f64; 8] = [
    0.189450610455069,
    0.182603415044924,
    0.169156519395003,
    0.149595988816577,
    0.124628971255534,
    0.095158511682493,
    0.062253523938648,
    0.027152459411754,
];
const LEGENDRE_WEIGHT_17: [f64; 9] = [
    0.179446470356206526,
    0.17656270536699265,
    0.168004102156450045,
    0.154045761076810288,
    0.13513636846852547,
    0.111883847193403971,
    0.085036148317179181,
    0.055459529373987201,
    0.024148302868547932,
];
const LEGENDRE_WEIGHT_18: [f64; 9] = [
    0.169142382963143592,
    0.16427648374583272,
    0.15468467512626524,
    0.140642914670650651,
    0.12255520671147846,
    0.100942044106287166,
    0.07642573025488906,
    0.049714548894969796,
    0.02161601352648331,
];
const LEGENDRE_WEIGHT_19: [f64; 10] = [
    0.161054449848783696,
    0.158968843393954348,
    0.152766042065859667,
    0.14260670217360661,
    0.128753962539336228,
    0.111566645547333995,
    0.091490021622449999,
    0.06904454273764123,
    0.0448142267656996,
    0.019461788229726477,
];
const LEGENDRE_WEIGHT_20: [f64; 10] = [
    0.15275338713072585,
    0.149172986472603747,
    0.142096109318382051,
    0.13168863844917663,
    0.118194531961518417,
    0.10193011981724044,
    0.083276741576704749,
    0.062672048334109064,
    0.040601429800386941,
    0.0176140071391521183,
];
const LEGENDRE_WEIGHT_21: [f64; 11] = [
    0.146081133649690427,
    0.144524403989970059,
    0.139887394791073155,
    0.13226893863333746,
    0.12183141605372853,
    0.108797299167148378,
    0.093444423456033862,
    0.0761001136283793,
    0.05713442542685721,
    0.036953789770852494,
    0.0160172282577743333,
];
const LEGENDRE_WEIGHT_22: [f64; 11] = [
    0.139251872855631993,
    0.136541498346015171,
    0.131173504787062371,
    0.123252376810512424,
    0.112932296080539218,
    0.100414144442880965,
    0.085941606217067727,
    0.06979646842452049,
    0.052293335152683286,
    0.033774901584814155,
    0.014627995298272201,
];
const LEGENDRE_WEIGHT_23: [f64; 12] = [
    0.13365457218610618,
    0.13246203940469662,
    0.12890572218808215,
    0.12304908430672953,
    0.114996640222411365,
    0.10489209146454141,
    0.09291576606003515,
    0.079281411776718955,
    0.064232421408525852,
    0.048037671731084669,
    0.030988005856979444,
    0.013411859487141772,
];
const LEGENDRE_WEIGHT_24: [f64; 12] = [
    0.127938195346752157,
    0.1258374563468283,
    0.121670472927803391,
    0.1155056680537256,
    0.10744427011596563,
    0.09761865210411389,
    0.08619016153195328,
    0.07334648141108031,
    0.0592985849154367807,
    0.04427743881741981,
    0.028531388628933663,
    0.0123412297999872,
];
const LEGENDRE_WEIGHT_25: [f64; 13] = [
    0.123176053726715451,
    0.122242442990310042,
    0.11945576353578477,
    0.114858259145711648,
    0.10851962447426365,
    0.10053594906705064,
    0.09102826198296365,
    0.08014070033500102,
    0.068038333812356917,
    0.054904695975835192,
    0.040939156701306313,
    0.026354986615032137,
    0.011393798501026288,
];
const LEGENDRE_WEIGHT_26: [f64; 13] = [
    0.118321415279262277,
    0.116660443485296582,
    0.11336181654631967,
    0.10847184052857659,
    0.10205916109442542,
    0.094213800355914148,
    0.085045894313485239,
    0.074684149765659746,
    0.063274046329574836,
    0.050975825297147812,
    0.03796238329436276,
    0.024417851092631909,
    0.01055137261734301,
];
const LEGENDRE_WEIGHT_27: [f64; 14] = [
    0.114220867378956989,
    0.11347634610896515,
    0.111252488356845193,
    0.10757828578853319,
    0.1025016378177458,
    0.09608872737002851,
    0.08842315854375695,
    0.079604867773057771,
    0.069748823766245593,
    0.0589835368598336,
    0.047449412520615063,
    0.03529705375741971,
    0.02268623159618062,
    0.00979899605129436,
];
const LEGENDRE_WEIGHT_28: [f64; 14] = [
    0.1100470130164752,
    0.108711192258294135,
    0.10605576592284642,
    0.10211296757806077,
    0.09693065799792992,
    0.0905717443930328409,
    0.083113417228901218,
    0.074646214234568779,
    0.065272923966999596,
    0.055107345675716745,
    0.04427293475900423,
    0.03290142778230438,
    0.02113211259277126,
    0.00912428259309452,
];
const LEGENDRE_WEIGHT_29: [f64; 15] = [
    0.106479381718314244,
    0.10587615509732094,
    0.104073310077729374,
    0.10109127375991497,
    0.096963834094408606,
    0.0917377571392587633,
    0.085472257366172528,
    0.0782383271357637838,
    0.07011793325505128,
    0.061203090657079139,
    0.05159482690249792,
    0.041402062518682836,
    0.030740492202093623,
    0.019732085056122706,
    0.00851690387874641,
];
const LEGENDRE_WEIGHT_30: [f64; 15] = [
    0.10285265289355884,
    0.1017623897484055,
    0.09959342058679527,
    0.09636873717464426,
    0.09212252223778613,
    0.08689978720108298,
    0.08075589522942022,
    0.073755974737705206,
    0.0659742298821805,
    0.057493156217619066,
    0.048402672830594053,
    0.03879919256962705,
    0.028784707883323369,
    0.018466468311090959,
    0.007968192496166606,
];

// =============================================================================
// Table for Kronrod Quadrature
// =============================================================================
const KRONROD_NODES_15: [f64; 8] = [
    0f64,
    0.207784955007898468,
    0.405845151377397167,
    0.58608723546769113,
    0.74153118559939444,
    0.864864423359769073,
    0.949107912342758525,
    0.991455371120812639,
];

const KRONROD_WEIGHTS_15: [f64; 8] = [
    0.20948214108472783,
    0.204432940075298892,
    0.19035057806478541,
    0.1690047266392679,
    0.140653259715525919,
    0.10479001032225018,
    0.06309209262997855,
    0.022935322010529225,
];

const KRONROD_NODES_21: [f64; 11] = [
    0f64,
    0.148874338981631211,
    0.294392862701460198,
    0.433395394129247191,
    0.562757134668604683,
    0.679409568299024406,
    0.780817726586416897,
    0.865063366688984511,
    0.930157491355708226,
    0.97390652851717172,
    0.995657163025808081,
];

const KRONROD_WEIGHTS_21: [f64; 11] = [
    0.149445554002916906,
    0.14773910490133849,
    0.142775938577060081,
    0.13470921731147333,
    0.12349197626206585,
    0.109387158802297642,
    0.09312545458369761,
    0.07503967481091995,
    0.054755896574351996,
    0.032558162307964727,
    0.011694638867371874,
];

const KRONROD_NODES_31: [f64; 16] = [
    0f64,
    0.101142066918717499,
    0.20119409399743452,
    0.29918000715316881,
    0.39415134707756337,
    0.485081863640239681,
    0.570972172608538848,
    0.650996741297416971,
    0.724417731360170047,
    0.790418501442465933,
    0.848206583410427216,
    0.897264532344081901,
    0.937273392400705904,
    0.967739075679139134,
    0.987992518020485429,
    0.99800229869339706,
];

const KRONROD_WEIGHTS_31: [f64; 16] = [
    0.101330007014791549,
    0.100769845523875595,
    0.09917359872179196,
    0.096642726983623679,
    0.093126598170825321,
    0.08856444305621177,
    0.083080502823133021,
    0.07684968075772038,
    0.06985412131872826,
    0.06200956780067064,
    0.053481524690928087,
    0.044589751324764877,
    0.0353463607913758462,
    0.02546084732671532,
    0.015007947329316123,
    0.005377479872923349,
];

const KRONROD_NODES_41: [f64; 21] = [
    0f64,
    0.076526521133497334,
    0.152605465240922676,
    0.227785851141645078,
    0.301627868114913004,
    0.373706088715419561,
    0.443593175238725103,
    0.510867001950827098,
    0.575140446819710315,
    0.636053680726515026,
    0.693237656334751385,
    0.746331906460150793,
    0.795041428837551198,
    0.839116971822218823,
    0.878276811252281976,
    0.912234428251325906,
    0.940822633831754754,
    0.963971927277913791,
    0.981507877450250259,
    0.993128599185094925,
    0.998859031588277664,
];

const KRONROD_WEIGHTS_41: [f64; 21] = [
    0.07660071191799966,
    0.076377867672080737,
    0.07570449768455667,
    0.074582875400499189,
    0.07303069033278667,
    0.07105442355344407,
    0.0686486729285216193,
    0.065834597133618422,
    0.062653237554781168,
    0.0591114008806395724,
    0.05519510534828599,
    0.050944573923728692,
    0.046434821867497675,
    0.04166887332797369,
    0.036600169758200798,
    0.0312873067770328,
    0.02588213360495116,
    0.020388373461266524,
    0.014626169256971253,
    0.008600269855642942,
    0.00307358371852053,
];

const KRONROD_NODES_51: [f64; 26] = [
    0f64,
    0.061544483005685079,
    0.122864692610710396,
    0.18371893942104889,
    0.243866883720988432,
    0.30308953893110783,
    0.361172305809387838,
    0.417885382193037749,
    0.47300273144571496,
    0.526325284334719183,
    0.577662930241222968,
    0.626810099010317413,
    0.673566368473468365,
    0.717766406813084388,
    0.759259263037357631,
    0.797873797998500059,
    0.833442628760834001,
    0.865847065293275595,
    0.894991997878275369,
    0.920747115281701562,
    0.942974571228974339,
    0.961614986425842512,
    0.976663921459517512,
    0.988035794534077248,
    0.995556969790498098,
    0.999262104992609834,
];

const KRONROD_WEIGHTS_51: [f64; 26] = [
    0.061580818067832935,
    0.061471189871425317,
    0.06112850971705305,
    0.060539455376045863,
    0.05972034032417406,
    0.058689680022394208,
    0.057437116361567833,
    0.055950811220412317,
    0.0542511298885454901,
    0.052362885806407476,
    0.050277679080715672,
    0.04798253713883671,
    0.04550291304992179,
    0.04287284502017005,
    0.040083825504032382,
    0.03711627148341554,
    0.03400213027432934,
    0.030792300167387489,
    0.027475317587851738,
    0.02400994560695322,
    0.02043537114588284,
    0.016847817709128298,
    0.013236229195571675,
    0.009473973386174152,
    0.005561932135356714,
    0.001987383892330316,
];

const KRONROD_NODES_61: [f64; 31] = [
    0f64,
    0.051471842555317696,
    0.10280693796673703,
    0.153869913608583547,
    0.204525116682309891,
    0.254636926167889846,
    0.30407320227362508,
    0.35270472553087811,
    0.40040125483039439,
    0.447033769538089177,
    0.492480467861778575,
    0.536624148142019899,
    0.579345235826361692,
    0.620526182989242861,
    0.66006106412662696,
    0.697850494793315797,
    0.733790062453226805,
    0.767777432104826195,
    0.799727835821839083,
    0.829565762382768397,
    0.857205233546061099,
    0.882560535792052682,
    0.905573307699907799,
    0.926200047429274326,
    0.944374444748559979,
    0.960021864968307512,
    0.973116322501126268,
    0.98366812327974721,
    0.991630996870404595,
    0.99689348407464954,
    0.999484410050490638,
];

const KRONROD_WEIGHTS_61: [f64; 31] = [
    0.05149472942945157,
    0.051426128537459026,
    0.051221547849258772,
    0.050881795898749606,
    0.0504059214027823468,
    0.049795683427074206,
    0.04905543455502978,
    0.048185861757087129,
    0.04718554656929915,
    0.046059238271006988,
    0.044814800133162663,
    0.043452539701356069,
    0.041969810215164246,
    0.04037453895153596,
    0.038678945624727593,
    0.036882364651821229,
    0.034979338028060024,
    0.03298144705748373,
    0.030907257562387762,
    0.028754048765041293,
    0.0265099548823331016,
    0.024191162078080601,
    0.02182803582160919,
    0.019414141193942381,
    0.016920889189053273,
    0.0143697295070458,
    0.011823015253496342,
    0.009273279659517763,
    0.006630703915931292,
    0.003890461127099884,
    0.001389013698677008,
];
