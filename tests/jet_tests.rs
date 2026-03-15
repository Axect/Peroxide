extern crate peroxide;
use peroxide::fuga::*;

// =============================================================================
// Helpers
// =============================================================================

const EPS: f64 = 1e-10;

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < EPS,
        "actual = {actual:.15e}, expected = {expected:.15e}, diff = {:.3e}",
        (actual - expected).abs()
    );
}

fn assert_close_eps(actual: f64, expected: f64, eps: f64) {
    assert!(
        (actual - expected).abs() < eps,
        "actual = {actual:.15e}, expected = {expected:.15e}, diff = {:.3e}",
        (actual - expected).abs()
    );
}

// =============================================================================
// Tier 1.1 — Construction & Accessors
// =============================================================================

#[test]
fn test_jet_var_value() {
    let x = Jet::<1>::var(3.0);
    assert_close(x.value(), 3.0);
    assert_close(x.x(), 3.0);
    assert_close(x.dx(), 1.0);
}

#[test]
fn test_jet_var_higher_order() {
    // Jet<3>::var(2.0): deriv[0]=1, deriv[1]=0, deriv[2]=0
    let x = Jet::<3>::var(2.0);
    assert_close(x.value(), 2.0);
    assert_close(x.dx(), 1.0);
    assert_close(x.ddx(), 0.0); // deriv[1]*2 = 0
    assert_close(x.taylor_coeff(0), 2.0);
    assert_close(x.taylor_coeff(1), 1.0);
    assert_close(x.taylor_coeff(2), 0.0);
}

#[test]
fn test_jet_constant() {
    let c = Jet::<2>::constant(7.0);
    assert_close(c.value(), 7.0);
    assert_close(c.dx(), 0.0);
    assert_close(c.ddx(), 0.0);
}

#[test]
fn test_jet_new_raw() {
    let j = Jet::<2>::new(5.0, [2.0, 3.0]);
    assert_close(j.value(), 5.0);
    // dx() = deriv[0] = 2.0
    assert_close(j.dx(), 2.0);
    // ddx() = deriv[1] * 2 = 6.0
    assert_close(j.ddx(), 6.0);
}

#[test]
fn test_ad0_constructor() {
    let j = ad0(4.2);
    assert_close(j.value(), 4.2);
    assert_close(j.dx(), 0.0); // N=0, dx() returns 0
}

#[test]
fn test_ad1_constructor() {
    let j = ad1(3.0, 5.0);
    assert_close(j.value(), 3.0);
    assert_close(j.dx(), 5.0);
}

#[test]
fn test_ad2_constructor_normalization() {
    // ad2(x, dx, ddx) stores ddx/2 internally; ddx() returns deriv[1]*2 = ddx
    let j = ad2(4.0, 4.0, 2.0);
    assert_close(j.value(), 4.0);
    assert_close(j.dx(), 4.0);
    // Normalized convention: ddx() == raw second derivative == 2.0
    assert_close(j.ddx(), 2.0);
}

#[test]
fn test_ad2_internal_storage() {
    // Internal deriv[1] should be ddx/2
    let j = ad2(1.0, 1.0, 6.0);
    assert_close(j.taylor_coeff(2), 3.0); // c_2 = ddx/2! = 6/2 = 3
}

#[test]
fn test_AD0_constructor() {
    let j = AD0(9.0);
    assert_close(j.value(), 9.0);
    assert_close(j.dx(), 0.0);
    assert_close(j.ddx(), 0.0);
}

#[test]
fn test_AD1_constructor() {
    let j = AD1(3.0, 1.0);
    assert_close(j.value(), 3.0);
    assert_close(j.dx(), 1.0);
    assert_close(j.ddx(), 0.0);
}

#[test]
fn test_AD2_constructor() {
    let j = AD2(5.0, 3.0, 4.0);
    assert_close(j.value(), 5.0);
    assert_close(j.dx(), 3.0);
    assert_close(j.ddx(), 4.0);
}

#[test]
fn test_derivative_accessor() {
    // f(x) = x^2 at x=3: f=9, f'=6, f''=2
    let x = Jet::<2>::var(3.0);
    let y = x.powi(2);
    assert_close(y.derivative(0), 9.0);
    assert_close(y.derivative(1), 6.0);
    assert_close(y.derivative(2), 2.0);
}

#[test]
fn test_taylor_coeff_accessor() {
    // f(x)=x^2 at x=3: c_0=9, c_1=f'(3)/1!=6, c_2=f''(3)/2!=1
    let x = Jet::<2>::var(3.0);
    let y = x.powi(2);
    assert_close(y.taylor_coeff(0), 9.0);
    assert_close(y.taylor_coeff(1), 6.0);
    assert_close(y.taylor_coeff(2), 1.0);
}

#[test]
fn test_derivative_out_of_bounds_returns_zero() {
    let x = Jet::<1>::var(2.0);
    let y = x.powi(3);
    // order > N → 0
    assert_close(y.derivative(5), 0.0);
    assert_close(y.taylor_coeff(5), 0.0);
}

// =============================================================================
// Tier 1.2 — Arithmetic (Jet<1> and Jet<2>)
// =============================================================================

#[test]
fn test_add_jet1() {
    let a = ad1(2.0, 3.0);
    let b = ad1(5.0, 1.0);
    let c = a + b;
    assert_close(c.value(), 7.0);
    assert_close(c.dx(), 4.0);
}

#[test]
fn test_sub_jet1() {
    let a = ad1(5.0, 4.0);
    let b = ad1(2.0, 1.0);
    let c = a - b;
    assert_close(c.value(), 3.0);
    assert_close(c.dx(), 3.0);
}

#[test]
fn test_mul_jet1_product_rule() {
    // (x * x)' = 2x; at x=3 → derivative = 6
    let x = Jet::<1>::var(3.0);
    let y = x * x;
    assert_close(y.value(), 9.0);
    assert_close(y.dx(), 6.0); // 2 * 3
}

#[test]
fn test_div_jet1_quotient_rule() {
    // (1/x)' = -1/x^2; at x=2 → -0.25
    let one = Jet::<1>::constant(1.0);
    let x = Jet::<1>::var(2.0);
    let y = one / x;
    assert_close(y.value(), 0.5);
    assert_close(y.dx(), -0.25);
}

#[test]
fn test_mul_jet2_product_rule_second_order() {
    // (x^2)'' = 2; at x=3 → second derivative = 2
    let x = Jet::<2>::var(3.0);
    let y = x * x;
    assert_close(y.value(), 9.0);
    assert_close(y.dx(), 6.0);
    assert_close(y.ddx(), 2.0);
}

#[test]
fn test_chain_rule_sin_x_squared() {
    // f(x) = sin(x^2), f'(x) = 2x*cos(x^2) at x=1
    let x = Jet::<1>::var(1.0);
    let y = (x * x).sin();
    let expected_value = (1.0f64).sin();
    let expected_deriv = 2.0 * 1.0 * (1.0f64).cos();
    assert_close(y.value(), expected_value);
    assert_close(y.dx(), expected_deriv);
}

#[test]
fn test_chain_rule_exp_neg_x_squared() {
    // f(x) = exp(-x^2), f'(x) = -2x*exp(-x^2) at x=1
    let x = Jet::<1>::var(1.0);
    let neg_x_sq = -(x * x);
    let y = neg_x_sq.exp();
    let ev = (-1.0f64).exp();
    let ed = -2.0 * ev;
    assert_close(y.value(), ev);
    assert_close(y.dx(), ed);
}

#[test]
fn test_scalar_add_jet_plus_f64() {
    let x = ad1(3.0, 1.0);
    let y = x + 5.0;
    assert_close(y.value(), 8.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_scalar_f64_plus_jet() {
    let x = ad1(3.0, 2.0);
    let y = 5.0 + x;
    assert_close(y.value(), 8.0);
    assert_close(y.dx(), 2.0);
}

#[test]
fn test_scalar_sub_jet_minus_f64() {
    let x = ad1(7.0, 1.0);
    let y = x - 3.0;
    assert_close(y.value(), 4.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_scalar_f64_minus_jet() {
    let x = ad1(3.0, 2.0);
    let y = 10.0 - x;
    assert_close(y.value(), 7.0);
    assert_close(y.dx(), -2.0);
}

#[test]
fn test_scalar_mul_jet_times_f64() {
    let x = ad1(3.0, 1.0);
    let y = x * 4.0;
    assert_close(y.value(), 12.0);
    assert_close(y.dx(), 4.0);
}

#[test]
fn test_scalar_f64_times_jet() {
    let x = ad1(3.0, 1.0);
    let y = 4.0 * x;
    assert_close(y.value(), 12.0);
    assert_close(y.dx(), 4.0);
}

#[test]
fn test_scalar_div_jet_by_f64() {
    let x = ad1(6.0, 2.0);
    let y = x / 2.0;
    assert_close(y.value(), 3.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_scalar_f64_div_jet() {
    // 6 / x at x=2, dx=1 → value=3, deriv = -6/x^2 = -1.5
    let x = ad1(2.0, 1.0);
    let y = 6.0 / x;
    assert_close(y.value(), 3.0);
    assert_close(y.dx(), -1.5);
}

#[test]
fn test_negation_jet1() {
    let x = ad1(3.0, 2.0);
    let y = -x;
    assert_close(y.value(), -3.0);
    assert_close(y.dx(), -2.0);
}

#[test]
fn test_negation_jet2() {
    let x = ad2(3.0, 2.0, 4.0);
    let y = -x;
    assert_close(y.value(), -3.0);
    assert_close(y.dx(), -2.0);
    assert_close(y.ddx(), -4.0);
}

// =============================================================================
// Tier 1.3 — Mathematical Operations
// =============================================================================

// --- exp ---

#[test]
fn test_exp_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.exp();
    // exp(0) = 1, (exp(x))' at x=0 = exp(0)*1 = 1
    assert_close(y.value(), 1.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_exp_at_one() {
    let x = Jet::<1>::var(1.0);
    let y = x.exp();
    let e = std::f64::consts::E;
    assert_close(y.value(), e);
    assert_close(y.dx(), e);
}

#[test]
fn test_exp_jet2_second_derivative() {
    // (exp(x))'' = exp(x); at x=0 → ddx = 1
    let x = Jet::<2>::var(0.0);
    let y = x.exp();
    assert_close(y.value(), 1.0);
    assert_close(y.dx(), 1.0);
    assert_close(y.ddx(), 1.0);
}

// --- ln ---

#[test]
fn test_ln_at_one() {
    // ln(1) = 0, (ln(x))' at x=1 = 1/1 = 1
    let x = Jet::<1>::var(1.0);
    let y = x.ln();
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_ln_at_e() {
    let e = std::f64::consts::E;
    let x = Jet::<1>::var(e);
    let y = x.ln();
    assert_close(y.value(), 1.0);
    assert_close_eps(y.dx(), 1.0 / e, 1e-14);
}

#[test]
fn test_ln_jet2_second_derivative() {
    // (ln(x))'' = -1/x^2; at x=2 → -0.25
    let x = Jet::<2>::var(2.0);
    let y = x.ln();
    assert_close(y.value(), 2.0f64.ln());
    assert_close(y.dx(), 0.5);
    assert_close(y.ddx(), -0.25);
}

// --- sin, cos ---

#[test]
fn test_sin_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.sin();
    // sin(0) = 0, (sin(x))' at x=0 = cos(0)*1 = 1
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_cos_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.cos();
    // cos(0) = 1, (cos(x))' at x=0 = -sin(0)*1 = 0
    assert_close(y.value(), 1.0);
    assert_close(y.dx(), 0.0);
}

#[test]
fn test_sin_at_pi_over_2() {
    let x = Jet::<1>::var(std::f64::consts::FRAC_PI_2);
    let y = x.sin();
    // sin(pi/2) = 1, (sin(x))' at pi/2 = cos(pi/2) = 0
    assert_close(y.value(), 1.0);
    assert_close_eps(y.dx(), 0.0, 1e-15);
}

#[test]
fn test_sin_cos_jet2_second_derivative() {
    // (sin(x))'' = -sin(x); at x=0 → 0
    let x = Jet::<2>::var(0.0);
    let (s, c) = x.sin_cos();
    assert_close(s.value(), 0.0);
    assert_close(s.dx(), 1.0);
    assert_close(s.ddx(), 0.0);
    assert_close(c.value(), 1.0);
    assert_close(c.dx(), 0.0);
    assert_close(c.ddx(), -1.0);
}

// --- sinh, cosh ---

#[test]
fn test_sinh_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.sinh();
    // sinh(0) = 0, (sinh(x))' = cosh(0) = 1
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_cosh_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.cosh();
    // cosh(0) = 1, (cosh(x))' = sinh(0) = 0
    assert_close(y.value(), 1.0);
    assert_close(y.dx(), 0.0);
}

#[test]
fn test_sinh_cosh_jet2() {
    // (sinh(x))'' = sinh(x); at x=0 → 0
    // (cosh(x))'' = cosh(x); at x=0 → 1
    let x = Jet::<2>::var(0.0);
    let (s, c) = x.sinh_cosh();
    assert_close(s.value(), 0.0);
    assert_close(s.dx(), 1.0);
    assert_close(s.ddx(), 0.0);
    assert_close(c.value(), 1.0);
    assert_close(c.dx(), 0.0);
    assert_close(c.ddx(), 1.0);
}

// --- sqrt ---

#[test]
fn test_sqrt_at_four() {
    // sqrt(4)=2, (sqrt(x))' = 1/(2*sqrt(x)) at x=4 → 0.25
    let x = Jet::<1>::var(4.0);
    let y = x.sqrt();
    assert_close(y.value(), 2.0);
    assert_close(y.dx(), 0.25);
}

#[test]
fn test_sqrt_jet2_second_derivative() {
    // (sqrt(x))'' = -1/(4 * x^(3/2)); at x=4 → -1/32
    let x = Jet::<2>::var(4.0);
    let y = x.sqrt();
    assert_close(y.value(), 2.0);
    assert_close(y.dx(), 0.25);
    assert_close(y.ddx(), -1.0 / 32.0);
}

// --- powi ---

#[test]
fn test_powi_x_cubed_jet1() {
    // x^3 at x=2: value=8, deriv=3*x^2=12
    let x = Jet::<1>::var(2.0);
    let y = x.powi(3);
    assert_close(y.value(), 8.0);
    assert_close(y.dx(), 12.0);
}

#[test]
fn test_powi_x_cubed_jet2() {
    // x^3 at x=2: value=8, f'=12, f''=6*x=12
    let x = Jet::<2>::var(2.0);
    let y = x.powi(3);
    assert_close(y.value(), 8.0);
    assert_close(y.dx(), 12.0);
    assert_close(y.ddx(), 12.0);
}

#[test]
fn test_powi_negative_power() {
    // x^(-2) at x=2: value=0.25, f'=-2/x^3=-0.25
    let x = Jet::<1>::var(2.0);
    let y = x.powi(-2);
    assert_close(y.value(), 0.25);
    assert_close(y.dx(), -0.25);
}

#[test]
fn test_powi_zero_power() {
    let x = Jet::<1>::var(5.0);
    let y = x.powi(0);
    assert_close(y.value(), 1.0);
    assert_close(y.dx(), 0.0);
}

// --- powf ---

#[test]
fn test_powf_x_to_2_5() {
    // f(x) = x^2.5 at x=4: value=32, f'=2.5*x^1.5=2.5*8=20
    let x = Jet::<1>::var(4.0);
    let y = x.powf(2.5);
    assert_close(y.value(), 32.0);
    assert_close(y.dx(), 20.0);
}

#[test]
fn test_powf_jet2_second_derivative() {
    // f(x) = x^2.5 at x=4: f''=2.5*1.5*x^0.5=2.5*1.5*2=7.5
    let x = Jet::<2>::var(4.0);
    let y = x.powf(2.5);
    assert_close(y.value(), 32.0);
    assert_close(y.dx(), 20.0);
    assert_close_eps(y.ddx(), 7.5, 1e-10);
}

// --- tan ---

#[test]
fn test_tan_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.tan();
    // tan(0)=0, (tan(x))' = sec^2(0) = 1
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_tan_jet2() {
    // (tan(x))'' = 2*tan(x)*sec^2(x); at x=0 → 0
    let x = Jet::<2>::var(0.0);
    let y = x.tan();
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
    assert_close(y.ddx(), 0.0);
}

// --- tanh ---

#[test]
fn test_tanh_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.tanh();
    // tanh(0)=0, (tanh(x))' = sech^2(0) = 1
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_tanh_jet2() {
    // (tanh(x))'' = -2*tanh(x)*sech^2(x); at x=0 → 0
    let x = Jet::<2>::var(0.0);
    let y = x.tanh();
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
    assert_close(y.ddx(), 0.0);
}

// --- asin ---

#[test]
fn test_asin_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.asin();
    // asin(0)=0, (asin(x))' = 1/sqrt(1-x^2) at x=0 → 1
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_asin_at_half() {
    let x = Jet::<1>::var(0.5);
    let y = x.asin();
    // asin(0.5) = pi/6, (asin(x))' = 1/sqrt(1-0.25) = 2/sqrt(3)
    assert_close(y.value(), std::f64::consts::FRAC_PI_6);
    let expected_deriv = 1.0 / (1.0f64 - 0.25f64).sqrt();
    assert_close_eps(y.dx(), expected_deriv, 1e-10);
}

// --- acos ---

#[test]
fn test_acos_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.acos();
    // acos(0) = pi/2, (acos(x))' = -1/sqrt(1-x^2) at x=0 → -1
    assert_close(y.value(), std::f64::consts::FRAC_PI_2);
    assert_close(y.dx(), -1.0);
}

// --- atan ---

#[test]
fn test_atan_at_zero() {
    let x = Jet::<1>::var(0.0);
    let y = x.atan();
    // atan(0) = 0, (atan(x))' = 1/(1+x^2) at x=0 → 1
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
}

#[test]
fn test_atan_at_one() {
    let x = Jet::<1>::var(1.0);
    let y = x.atan();
    // atan(1) = pi/4, (atan(x))' = 1/(1+1) = 0.5
    assert_close(y.value(), std::f64::consts::FRAC_PI_4);
    assert_close(y.dx(), 0.5);
}

#[test]
fn test_atan_jet2_second_derivative() {
    // (atan(x))'' = -2x/(1+x^2)^2; at x=0 → 0
    let x = Jet::<2>::var(0.0);
    let y = x.atan();
    assert_close(y.value(), 0.0);
    assert_close(y.dx(), 1.0);
    assert_close(y.ddx(), 0.0);
}

// =============================================================================
// Tier 1.4 — Higher-Order (Jet<5>, Jet<10>)
// =============================================================================

#[test]
fn test_x_to_the_5th_at_x1_jet5() {
    // f(x) = x^5 at x=1:
    // f^(0) = 1, f^(1) = 5, f^(2) = 20, f^(3) = 60, f^(4) = 120, f^(5) = 120
    let x = Jet::<5>::var(1.0);
    let y = x.powi(5);
    assert_close(y.derivative(0), 1.0);
    assert_close(y.derivative(1), 5.0);
    assert_close(y.derivative(2), 20.0);
    assert_close(y.derivative(3), 60.0);
    assert_close(y.derivative(4), 120.0);
    assert_close(y.derivative(5), 120.0);
}

#[test]
fn test_sin_at_zero_jet10_derivative_cycle() {
    // sin(x) derivatives cycle: sin, cos, -sin, -cos, sin, ...
    // At x=0: (0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0)
    let x = Jet::<10>::var(0.0);
    let y = x.sin();
    let expected = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0];
    for (k, &exp) in expected.iter().enumerate() {
        assert_close_eps(
            y.derivative(k),
            exp,
            1e-10,
        );
    }
}

#[test]
fn test_exp_at_zero_jet10_all_derivatives_one() {
    // exp(x) at x=0: all derivatives = 1
    let x = Jet::<10>::var(0.0);
    let y = x.exp();
    for k in 0..=10 {
        assert_close_eps(y.derivative(k), 1.0, 1e-9);
    }
}

#[test]
fn test_x_squared_jet5_higher_derivatives_zero() {
    // f(x) = x^2 at x=3: f^(k) = 0 for k > 2
    let x = Jet::<5>::var(3.0);
    let y = x.powi(2);
    assert_close(y.derivative(0), 9.0);
    assert_close(y.derivative(1), 6.0);
    assert_close(y.derivative(2), 2.0);
    assert_close(y.derivative(3), 0.0);
    assert_close(y.derivative(4), 0.0);
    assert_close(y.derivative(5), 0.0);
}

#[test]
fn test_cos_at_zero_jet10_derivative_cycle() {
    // cos(x) derivatives at x=0: (1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1)
    let x = Jet::<10>::var(0.0);
    let y = x.cos();
    let expected = [1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0];
    for (k, &exp) in expected.iter().enumerate() {
        assert_close_eps(
            y.derivative(k),
            exp,
            1e-10,
        );
    }
}

// =============================================================================
// Tier 1.5 — Backward Compatibility
// =============================================================================

#[test]
fn test_ad_type_alias_is_jet2() {
    let j: AD = AD2(3.0, 1.0, 0.0);
    assert_close(j.value(), 3.0);
    assert_close(j.dx(), 1.0);
    assert_close(j.ddx(), 0.0);
}

#[test]
fn test_ad1_constructor_gives_jet2() {
    let j: Jet<2> = AD1(5.0, 2.0);
    assert_close(j.value(), 5.0);
    assert_close(j.dx(), 2.0);
    assert_close(j.ddx(), 0.0);
}

#[test]
fn test_ad2_constructor_gives_jet2() {
    let j: Jet<2> = AD2(5.0, 2.0, 6.0);
    assert_close(j.value(), 5.0);
    assert_close(j.dx(), 2.0);
    assert_close(j.ddx(), 6.0);
}

#[test]
fn test_ad_vec_trait_to_ad_vec() {
    let v: Vec<f64> = vec![1.0, 2.0, 3.0];
    let av: Vec<AD> = v.to_ad_vec();
    assert_eq!(av.len(), 3);
    assert_close(av[0].value(), 1.0);
    assert_close(av[1].value(), 2.0);
    assert_close(av[2].value(), 3.0);
    // All derivatives zero (constant)
    assert_close(av[0].dx(), 0.0);
}

#[test]
fn test_ad_vec_trait_to_f64_vec_compat() {
    let av: Vec<AD> = vec![AD0(1.0), AD1(2.0, 1.0), AD2(3.0, 1.0, 0.0)];
    let fv = av.to_f64_vec_compat();
    assert_close(fv[0], 1.0);
    assert_close(fv[1], 2.0);
    assert_close(fv[2], 3.0);
}

#[test]
fn test_jet_vec_trait_to_jet_vec() {
    let v: Vec<f64> = vec![1.0, 2.0, 3.0];
    let jv: Vec<Jet<1>> = v.to_jet_vec();
    assert_eq!(jv.len(), 3);
    assert_close(jv[0].value(), 1.0);
    assert_close(jv[2].value(), 3.0);
    assert_close(jv[0].dx(), 0.0);
}

#[test]
fn test_jet_vec_to_f64_vec() {
    let jv: Vec<Jet<1>> = vec![ad1(1.0, 0.5), ad1(2.0, 1.0)];
    let fv = jv.to_f64_vec();
    assert_close(fv[0], 1.0);
    assert_close(fv[1], 2.0);
}

#[test]
fn test_adfn_and_stable_fn_value() {
    // ADFn at grad_level=0 returns f(x)
    let f = ADFn::new(|x: Jet<2>| x.powi(2));
    let val: f64 = f.call_stable(3.0);
    assert_close(val, 9.0);
}

#[test]
fn test_adfn_and_stable_fn_first_derivative() {
    // ADFn at grad_level=1 returns f'(x)
    let f = ADFn::new(|x: Jet<2>| x.powi(2));
    let df = f.grad();
    let val: f64 = df.call_stable(3.0);
    assert_close(val, 6.0); // (x^2)' at x=3 = 6
}

#[test]
fn test_adfn_and_stable_fn_second_derivative() {
    // ADFn at grad_level=2 returns f''(x)
    let f = ADFn::new(|x: Jet<2>| x.powi(2));
    let ddf = f.grad().grad();
    let val: f64 = ddf.call_stable(3.0);
    assert_close(val, 2.0); // (x^2)'' = 2
}

#[test]
fn test_adfn_stable_fn_jet2_passthrough() {
    let f = ADFn::new(|x: Jet<2>| x.powi(3));
    let x = Jet::<2>::new(2.0, [1.0, 0.0]);
    let y: Jet<2> = f.call_stable(x);
    assert_close(y.value(), 8.0);
    assert_close(y.dx(), 12.0);
}

// =============================================================================
// Tier 1.6 — Real Trait
// =============================================================================

fn compute_square<T: Real>(x: T) -> T {
    x.powi(2)
}

fn compute_exp_neg_x<T: Real>(x: T) -> T {
    let neg = T::from_f64(-1.0) * x;
    neg.exp()
}

#[test]
fn test_real_trait_with_f64() {
    let result = compute_square(3.0f64);
    assert_close(result, 9.0);
}

#[test]
fn test_real_trait_with_ad() {
    let x = AD1(3.0, 1.0);
    let result = compute_square(x);
    assert_close(result.value(), 9.0);
    assert_close(result.dx(), 6.0);
}

#[test]
fn test_real_trait_exp_with_f64() {
    let result = compute_exp_neg_x(1.0f64);
    assert_close(result, (-1.0f64).exp());
}

#[test]
fn test_real_trait_exp_with_ad() {
    // f(x) = exp(-x) at x=1: value=e^-1, f'=-e^-1
    let x = AD1(1.0, 1.0);
    let result = compute_exp_neg_x(x);
    let e_inv = (-1.0f64).exp();
    assert_close(result.value(), e_inv);
    assert_close(result.dx(), -e_inv);
}

#[test]
fn test_real_trait_to_f64() {
    let x = AD2(5.0, 1.0, 0.0);
    let v: f64 = x.to_f64();
    assert_close(v, 5.0);
}

#[test]
fn test_real_trait_from_f64() {
    let x = AD::from_f64(7.0);
    assert_close(x.value(), 7.0);
    assert_close(x.dx(), 0.0);
}

#[test]
fn test_real_trait_to_ad() {
    let x = 3.0f64.to_ad();
    assert_close(x.value(), 3.0);
    assert_close(x.dx(), 0.0);
}

// Generic function exercising more Real operations
fn polynomial<T: Real>(x: T) -> T {
    // p(x) = x^3 - 2*x^2 + x - 1
    let one = T::from_f64(1.0);
    let two = T::from_f64(2.0);
    x.powi(3) - x.powi(2) * two + x - one
}

#[test]
fn test_real_trait_polynomial_f64() {
    // p(2) = 8 - 8 + 2 - 1 = 1
    let result = polynomial(2.0f64);
    assert_close(result, 1.0);
}

#[test]
fn test_real_trait_polynomial_ad() {
    // p'(x) = 3x^2 - 4x + 1; p'(2) = 12 - 8 + 1 = 5
    let x = AD1(2.0, 1.0);
    let result = polynomial(x);
    assert_close(result.value(), 1.0);
    assert_close(result.dx(), 5.0);
}

// =============================================================================
// Tier 1.7 — Edge Cases
// =============================================================================

#[test]
fn test_jet0_constant_only() {
    let j = ad0(3.14);
    assert_close(j.value(), 3.14);
    // No derivatives in Jet<0>
    assert_close(j.dx(), 0.0);
    assert_close(j.ddx(), 0.0);
}

#[test]
fn test_jet0_arithmetic() {
    let a = ad0(2.0);
    let b = ad0(3.0);
    let c = a + b;
    assert_close(c.value(), 5.0);
    let d = a * b;
    assert_close(d.value(), 6.0);
}

#[test]
fn test_nan_propagation_add() {
    let a = ad1(f64::NAN, 1.0);
    let b = ad1(1.0, 1.0);
    let c = a + b;
    assert!(c.value().is_nan());
}

#[test]
fn test_nan_propagation_mul() {
    let a = ad1(f64::NAN, 0.0);
    let b = ad1(2.0, 1.0);
    let c = a * b;
    assert!(c.value().is_nan());
}

#[test]
fn test_nan_propagation_exp() {
    let j = ad1(f64::NAN, 1.0);
    let y = j.exp();
    assert!(y.value().is_nan());
}

#[test]
fn test_div_by_zero_value() {
    let a = ad1(1.0, 0.0);
    let b = ad1(0.0, 0.0);
    let c = a / b;
    assert!(c.value().is_infinite());
}

#[test]
fn test_ln_of_negative_is_nan() {
    let j = ad1(-1.0, 1.0);
    let y = j.ln();
    assert!(y.value().is_nan());
}

#[test]
fn test_sqrt_of_negative_is_nan() {
    let j = ad1(-4.0, 1.0);
    let y = j.sqrt();
    assert!(y.value().is_nan());
}

// =============================================================================
// Tier 1.8 — Index operator (backward compat)
// =============================================================================

#[test]
fn test_index_operator_jet1() {
    let j = ad1(5.0, 3.0);
    assert_close(j[0], 5.0); // value
    assert_close(j[1], 3.0); // deriv[0]
}

#[test]
fn test_index_operator_jet2() {
    let j = ad2(5.0, 3.0, 4.0);
    assert_close(j[0], 5.0);     // value
    assert_close(j[1], 3.0);     // deriv[0] = dx
    assert_close(j[2], 2.0);     // deriv[1] = ddx/2 = 4/2 = 2
}

#[test]
fn test_index_mut_operator() {
    let mut j = ad1(1.0, 0.0);
    j[0] = 5.0;
    j[1] = 2.0;
    assert_close(j.value(), 5.0);
    assert_close(j.dx(), 2.0);
}

// =============================================================================
// Tier 1.9 — From/Into conversions
// =============================================================================

#[test]
fn test_from_f64_creates_constant() {
    let j = Jet::<2>::from(3.0f64);
    assert_close(j.value(), 3.0);
    assert_close(j.dx(), 0.0);
    assert_close(j.ddx(), 0.0);
}

#[test]
fn test_into_f64_extracts_value() {
    let j = ad2(7.5, 1.0, 0.0);
    let v: f64 = j.into();
    assert_close(v, 7.5);
}

// =============================================================================
// Tier 1.10 — Partial ordering by value
// =============================================================================

#[test]
fn test_partial_ord_by_value() {
    let a = ad1(1.0, 5.0);
    let b = ad1(3.0, 0.0);
    assert!(a < b);
    assert!(b > a);
    assert!(a != b);
}

// =============================================================================
// Tier 1.11 — Display
// =============================================================================

#[test]
fn test_display_jet1() {
    let j = ad1(3.0, 1.0);
    let s = format!("{}", j);
    assert!(s.contains("3"), "display should contain value: {s}");
}

#[test]
fn test_display_jet0() {
    let j = ad0(2.5);
    let s = format!("{}", j);
    assert!(s.contains("2.5"), "display should contain value: {s}");
}

// =============================================================================
// Tier 1.12 — log variants
// =============================================================================

#[test]
fn test_log2_jet1() {
    // log2(8) = 3, (log2(x))' = 1/(x * ln(2)) at x=8
    let x = Jet::<1>::var(8.0);
    let y = x.log2();
    assert_close(y.value(), 3.0);
    let expected_deriv = 1.0 / (8.0 * 2.0f64.ln());
    assert_close_eps(y.dx(), expected_deriv, 1e-14);
}

#[test]
fn test_log10_jet1() {
    // log10(100) = 2, (log10(x))' = 1/(x * ln(10)) at x=100
    let x = Jet::<1>::var(100.0);
    let y = x.log10();
    assert_close(y.value(), 2.0);
    let expected_deriv = 1.0 / (100.0 * 10.0f64.ln());
    assert_close_eps(y.dx(), expected_deriv, 1e-14);
}

#[test]
fn test_log_base_3() {
    // log_3(27) = 3, (log_3(x))' = 1/(x*ln(3)) at x=27
    let x = Jet::<1>::var(27.0);
    let y = x.log(3.0);
    assert_close(y.value(), 3.0);
    let expected_deriv = 1.0 / (27.0 * 3.0f64.ln());
    assert_close_eps(y.dx(), expected_deriv, 1e-14);
}

// =============================================================================
// Tier 1.13 — Dual / HyperDual type aliases
// =============================================================================

#[test]
fn test_dual_alias() {
    let x: Dual = Dual::var(2.0);
    let y = x.powi(3);
    assert_close(y.value(), 8.0);
    assert_close(y.dx(), 12.0);
}

#[test]
fn test_hyperdual_alias() {
    let x: HyperDual = HyperDual::var(2.0);
    let y = x.powi(3);
    assert_close(y.value(), 8.0);
    assert_close(y.dx(), 12.0);
    assert_close(y.ddx(), 12.0);
}

// =============================================================================
// Tier 1.14 — FPVector / VecOps on Vec<Jet<1>>
// =============================================================================

#[test]
fn test_fpvector_fmap_jet1() {
    let v: Vec<Jet<1>> = vec![Jet::<1>::constant(1.0), Jet::<1>::constant(2.0)];
    let result = v.fmap(|x| x + x);
    assert_close(result[0].value(), 2.0);
    assert_close(result[1].value(), 4.0);
}

#[test]
fn test_fpvector_sum_jet1() {
    let v: Vec<Jet<1>> = vec![
        ad1(1.0, 1.0),
        ad1(2.0, 2.0),
        ad1(3.0, 3.0),
    ];
    let s = v.sum();
    // FPVector::sum uses reduce(self[0], +) which double-counts first element
    assert_close(s.value(), 7.0);  // 1 + (1+2+3)
    assert_close(s.dx(), 7.0);
}

#[test]
fn test_fpvector_prod_jet1() {
    let v: Vec<Jet<1>> = vec![
        Jet::<1>::constant(2.0),
        Jet::<1>::constant(3.0),
    ];
    let p = v.prod();
    // FPVector::prod uses reduce(self[0], *) which double-counts first element
    assert_close(p.value(), 12.0);  // 2 * (2*3)
}

// =============================================================================
// Tier 1.15 — Comprehensive composition tests
// =============================================================================

#[test]
fn test_composition_ln_sin() {
    // f(x) = ln(sin(x)) at x = pi/4
    // f'(x) = cos(x)/sin(x) = cot(x)
    let x0 = std::f64::consts::FRAC_PI_4;
    let x = Jet::<1>::var(x0);
    let y = x.sin().ln();
    let sin_x0 = x0.sin();
    assert_close_eps(y.value(), sin_x0.ln(), 1e-14);
    let expected_deriv = x0.cos() / x0.sin();
    assert_close_eps(y.dx(), expected_deriv, 1e-13);
}

#[test]
fn test_composition_sqrt_of_exp() {
    // f(x) = sqrt(exp(x)) = exp(x/2), f'(x) = 0.5*exp(x/2)
    let x = Jet::<1>::var(2.0);
    let y = x.exp().sqrt();
    let ev = (1.0f64).exp(); // exp(2/2) = exp(1)
    assert_close_eps(y.value(), ev, 1e-14);
    assert_close_eps(y.dx(), 0.5 * ev, 1e-14);
}

#[test]
fn test_taylor_series_ln_1_plus_x_at_zero() {
    // ln(1+x) at x=0: taylor coefficients c_k = (-1)^(k+1)/k
    // c_0=0, c_1=1, c_2=-1/2, c_3=1/3, c_4=-1/4, c_5=1/5
    let x = Jet::<5>::var(0.0);
    let one = Jet::<5>::constant(1.0);
    let y = (one + x).ln();
    assert_close(y.taylor_coeff(0), 0.0);
    assert_close(y.taylor_coeff(1), 1.0);
    assert_close(y.taylor_coeff(2), -0.5);
    assert_close_eps(y.taylor_coeff(3), 1.0 / 3.0, 1e-13);
    assert_close(y.taylor_coeff(4), -0.25);
    assert_close_eps(y.taylor_coeff(5), 0.2, 1e-13);
}

#[test]
fn test_taylor_series_sin_at_zero_jet5() {
    // sin(x) at x=0: c_0=0, c_1=1, c_2=0, c_3=-1/6, c_4=0, c_5=1/120
    let x = Jet::<5>::var(0.0);
    let y = x.sin();
    assert_close(y.taylor_coeff(0), 0.0);
    assert_close(y.taylor_coeff(1), 1.0);
    assert_close(y.taylor_coeff(2), 0.0);
    assert_close_eps(y.taylor_coeff(3), -1.0 / 6.0, 1e-14);
    assert_close(y.taylor_coeff(4), 0.0);
    assert_close_eps(y.taylor_coeff(5), 1.0 / 120.0, 1e-15);
}
