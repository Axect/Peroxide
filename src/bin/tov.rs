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
    let init_val = c!(0, 0, p_c);
    let records = solve_with_condition(tov_rhs, init_val, (0f64, 16f64), 1e-3, GL4(1e-15), |xs| xs[2] >= 0f64);
    records.print();
    records.write_with_header("example_data/tov_gl4.csv", vec!["r", "m", "phi", "p"]).expect("Can't write file");
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
pub fn tov_rhs(r: Dual, rs: Vec<Dual>) -> Vec<Dual> {
    let m_old = rs[0];
    let _phi_old = rs[1];
    let p_old = rs[2];

    let rho_0 = (p_old / K).powf(1. / GAMMAF);
    let rho_old = rho_0 + p_old / (GAMMAF - 1f64);

    let dm = 4f64 * PI * r.pow(2) * rho_old;

    let dphi = if r.value() != 0f64 { (m_old + 4f64 * PI * r.pow(3) * p_old) / (r.pow(2) * (1f64 - 2f64 * m_old / r)) } else { dual(0, 0) };
    let dp = if r.value() != 0f64 {- (rho_old + p_old) * (m_old + 4f64 * PI * r.pow(3) * p_old) / (r.pow(2) * (1f64 - 2f64 * m_old / r)) } else { dual(0, 0) };

    vec![dm, dphi, dp]
}