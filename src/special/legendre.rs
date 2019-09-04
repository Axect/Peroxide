use structure::polynomial::{poly, Polynomial};
use zip_with;

/// Legendre Polynomial
///
/// # Description
/// : Generate `n`-th order of Legendre polynomial
pub fn legendre_polynomial(n: usize) -> Polynomial {
    match n {
        0 => poly(vec![1f64]),       // 1
        1 => poly(vec![1f64, 0f64]), // x
        2 => poly(vec![1.5, 0f64, -0.5]),
        3 => poly(vec![2.5, 0f64, -1.5, 0f64]),
        _ => {
            let k = n - 1;
            let k_f64 = k as f64;
            ((2f64 * k_f64 + 1f64) * poly(vec![1f64, 0f64]) * legendre_polynomial(k)
                - k_f64 * legendre_polynomial(k - 1))
                / (k_f64 + 1f64)
        }
    }
}

pub fn unit_gauss_legendre_quadrature<F>(f: F, n: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    let (a, x) = gauss_legendre_table(n);
    let mut s = 0f64;
    for i in 0..a.len() {
        s += a[i] * f(x[i]);
    }
    s
}

/// Gauss Legendre Quadrature
///
/// # Type
/// * `f, n, (a,b) -> f64`
///     * `f`: Numerical function (`Fn(f64) -> f64`)
///     * `n`: Order of Legendre polynomial (up to 16)
///     * `(a,b)`: Interval of integration
///
/// # Reference
/// A. N. Lowan et al. (1942)
pub fn gauss_legendre_quadrature<F>(f: F, n: usize, (a, b): (f64, f64)) -> f64
where
    F: Fn(f64) -> f64,
{
    (b - a) / 2f64 * unit_gauss_legendre_quadrature(|x| f(x * (b - a) / 2f64 + (a + b) / 2f64), n)
}

pub fn gauss_legendre_table(n: usize) -> (Vec<f64>, Vec<f64>) {
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
        _ => unreachable!()
    }
    (result_weight, result_root)
}

// =============================================================================
// Table for Gauss-Legendre Quadrature (ref. A. N. Lowan et al. (1942))
// =============================================================================
pub const LEGENDRE_ROOT_2: [f64; 1] = [0.577350269189626];
pub const LEGENDRE_ROOT_3: [f64; 2] = [0f64, 0.774596669241483];
pub const LEGENDRE_ROOT_4: [f64; 2] = [0.339981043584856, 0.861136311594053];
pub const LEGENDRE_ROOT_5: [f64; 3] = [0f64, 0.538469310105683, 0.906179845938664];
pub const LEGENDRE_ROOT_6: [f64; 3] = [0.238619186083197, 0.661209386466265, 0.932469514203152];
pub const LEGENDRE_ROOT_7: [f64; 4] = [
    0f64,
    0.405845151377397,
    0.741531185599394,
    0.949107912342759,
];
pub const LEGENDRE_ROOT_8: [f64; 4] = [
    0.183434642495650,
    0.525532409916329,
    0.796666477413627,
    0.960289856497536,
];
pub const LEGENDRE_ROOT_9: [f64; 5] = [
    0f64,
    0.324253423403809,
    0.613371432700590,
    0.836031107326636,
    0.968160239507626,
];
pub const LEGENDRE_ROOT_10: [f64; 5] = [
    0.148874338981631,
    0.433395394129247,
    0.679409568299024,
    0.865063366688985,
    0.973906528517172,
];
pub const LEGENDRE_ROOT_11: [f64; 6] = [
    0f64,
    0.269543155952345,
    0.519096129110681,
    0.730152005574049,
    0.887062599768095,
    0.978228658146057,
];
pub const LEGENDRE_ROOT_12: [f64; 6] = [
    0.125333408511469,
    0.367831498918180,
    0.587317954286617,
    0.769902674194305,
    0.904117256370475,
    0.981560634246719,
];
pub const LEGENDRE_ROOT_13: [f64; 7] = [
    0f64,
    0.230458315955135,
    0.448492751036447,
    0.642349339440340,
    0.801578090733310,
    0.917598399222978,
    0.984183054718588,
];
pub const LEGENDRE_ROOT_14: [f64; 7] = [
    0.108054948707344,
    0.319112368927890,
    0.515248636358154,
    0.687292904811685,
    0.827201315069765,
    0.928434883663574,
    0.986283808696812,
];
pub const LEGENDRE_ROOT_15: [f64; 8] = [
    0f64,
    0.201194093997435,
    0.394151347077563,
    0.570972172608539,
    0.724417731360170,
    0.848206583410427,
    0.937273392400706,
    0.987992518020485,
];
pub const LEGENDRE_ROOT_16: [f64; 8] = [
    0.095012509837637,
    0.281603550779259,
    0.458016777657227,
    0.617876244402644,
    0.755404408355003,
    0.865631202387832,
    0.944575023073233,
    0.989400934991650,
];
pub const LEGENDRE_WEIGHT_2: [f64; 1] = [1f64];
pub const LEGENDRE_WEIGHT_3: [f64; 2] = [0.888888888888889, 0.555555555555556];
pub const LEGENDRE_WEIGHT_4: [f64; 2] = [0.652145154862546, 0.347854845137454];
pub const LEGENDRE_WEIGHT_5: [f64; 3] = [0.568888888888889, 0.478628670499366, 0.236926885056189];
pub const LEGENDRE_WEIGHT_6: [f64; 3] = [0.467913934572691, 0.360761573048139, 0.171324492379170];
pub const LEGENDRE_WEIGHT_7: [f64; 4] = [
    0.417959183673469,
    0.381830050505119,
    0.279705391489277,
    0.129484966168870,
];
pub const LEGENDRE_WEIGHT_8: [f64; 4] = [
    0.362683783378362,
    0.313706645877887,
    0.222381034453374,
    0.101228536290376,
];
pub const LEGENDRE_WEIGHT_9: [f64; 5] = [
    0.330239355001260,
    0.312347077040003,
    0.260610696402935,
    0.180648160694857,
    0.081274388361574,
];
pub const LEGENDRE_WEIGHT_10: [f64; 5] = [
    0.295524224714753,
    0.269266719309996,
    0.219086362515982,
    0.149451349150581,
    0.066671344308688,
];
pub const LEGENDRE_WEIGHT_11: [f64; 6] = [
    0.272925086777901,
    0.262804544510247,
    0.233193764591990,
    0.186290210927734,
    0.125580369464905,
    0.055668567116174,
];
pub const LEGENDRE_WEIGHT_12: [f64; 6] = [
    0.249147045813403,
    0.233492536538355,
    0.203167426723066,
    0.160078328543346,
    0.106939325995318,
    0.047175336386512,
];
pub const LEGENDRE_WEIGHT_13: [f64; 7] = [
    0.232551553230874,
    0.226283180262897,
    0.207816047536889,
    0.178145980761946,
    0.138873510219787,
    0.092121499837728,
    0.040484004765316,
];
pub const LEGENDRE_WEIGHT_14: [f64; 7] = [
    0.215263853463158,
    0.205198463721290,
    0.185538397477938,
    0.157203167158194,
    0.121518570687903,
    0.080158087159760,
    0.035119460331752,
];
pub const LEGENDRE_WEIGHT_15: [f64; 8] = [
    0.202578241925561,
    0.198431485327111,
    0.186161000015562,
    0.166269205816994,
    0.139570677926154,
    0.107159220467172,
    0.070366047488108,
    0.030753241996117,
];
pub const LEGENDRE_WEIGHT_16: [f64; 8] = [
    0.189450610455069,
    0.182603415044924,
    0.169156519395003,
    0.149595988816577,
    0.124628971255534,
    0.095158511682493,
    0.062253523938648,
    0.027152459411754,
];
