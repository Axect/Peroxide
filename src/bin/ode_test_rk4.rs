extern crate peroxide;
use peroxide::*;

use std::f64::consts::PI;

// ===========================================================
// Declare Constants
// ===========================================================
pub const K: f64 = 100.;
pub const GAMMA: usize = 2;
pub const GAMMAF: f64 = 2.;
pub const RHO0C: f64 = 0.4;

/// Correspond to a density of `2.2e14 g/cm^3` for the core limit
pub const RHO0L: f64 = 0.1324;

/// Correspond to ad density of `4.3e11 g/cm^3` for neutron drip
pub const RHO0D: f64 = 2.573e-4;

pub fn main() {
    let p_c = K*RHO0C.powf(GAMMAF);
    let init_val = c!(1e-15, 0, p_c);
    let init_state = State::new(0f64, init_val, vec![0f64; 3]);
    let mut ode_solver = ExplicitODE::new(tov);
    ode_solver
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state)
        .set_step_size(0.001)
        .set_times(10000);

    let results = ode_solver.integrate();

    let mut swt = SimpleWriter::new();
    swt
        .insert_matrix(results)
        .set_path("example_data/tov_rk4.pickle")
        .write_pickle();

}


/// Tolman-Oppenheimer-Volkoff Equations for Spherically Symmetric Equilibrium Stars
///
/// # Equation
/// ```latex
/// dm/dr = 4πr^2ρ
/// dP/dr = -ρm/r^2 (1 + P/ρ) (1 + 4πPr^3/m) (1 - 2m/r)^{-1}
/// dΦ/dr = -1/ρ dP/dr (1 + P/ρ)^{-1}
/// P = Kρ_0^Γ
/// Γ = 1 + 1/n
/// ρ = ρ_0(1 + ε)
/// m = 0, P = P_c at r = 0
/// ```
pub fn tov(st: &mut State<f64>) {
    let r = st.param;
    let rs = &st.state;
    let drs = &mut st.deriv;

    let m_old = rs[0];
    let _phi_old = rs[1];
    let p_old = rs[2];

    let rho_0 = (p_old / K).powf(1. / GAMMAF);
    let rho_old = rho_0 + p_old / (GAMMAF - 1f64);

    let dm = 4f64 * PI * r.powi(2) * rho_old;
    let dphi = if r != 0f64 { (m_old + 4f64 * PI * r.powi(3) * p_old) / (r.powi(2) * (1f64 - 2f64 * m_old / r)) } else { 0.0 };
    let dp = if r != 0f64 {- (rho_old + p_old) * (m_old + 4f64 * PI * r.powi(3) * p_old) / (r.powi(2) * (1f64 - 2f64 * m_old / r)) } else { 0.0 };

    drs[0] = dm;
    drs[1] = dphi;
    drs[2] = dp;
}