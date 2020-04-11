//! To find Eigenvalue & Eigenvector
//!
//! * Reference : Press, William H., and William T. Vetterling. *Numerical Recipes.* Cambridge: Cambridge Univ. Press, 2007.

#[cfg(feature = "O3")]
use lapack::{dsytrd, dorgtr, dsteqr};

use ::{Matrix, eye_shape};
use std::f64::EPSILON;
pub use self::EigenMethod::*;

#[derive(Debug, Copy, Clone)]
pub enum EigenMethod {
    Jacobi
}

#[derive(Debug, Clone)]
pub struct Eigen {
    pub eigenvalue: Vec<f64>,
    pub eigenvector: Matrix,
}

impl Eigen {
    pub fn extract(self) -> (Vec<f64>, Matrix) {
        (self.eigenvalue, self.eigenvector)
    }
}

pub fn eigen(m: &Matrix, em: EigenMethod) -> Eigen {
    match em {
        Jacobi => {
            let mat = m.clone();
            let mut j = jacobi(mat);
            j.iter();
            j.extract()
        }
    }
}

// =============================================================================
// Jacobi Method
// =============================================================================
/// To do Jacobi method
///
/// * Reference : Press, William H., and William T. Vetterling. *Numerical Recipes.* Cambridge: Cambridge Univ. Press, 2007.
#[derive(Debug)]
pub struct JacobiTemp {
    pub n: usize,
    pub a: Matrix,
    pub v: Matrix,
    pub d: Vec<f64>,
    pub n_rot: usize,
}

pub fn jacobi(m: Matrix) -> JacobiTemp {
    let n = m.row;
    let v = eye_shape(n, m.shape);
    let d = m.diag();
    let a = m;

    JacobiTemp {
        n,
        a,
        v,
        d,
        n_rot: 0
    }
}

impl JacobiTemp {
    pub fn new(m: Matrix) -> Self {
        jacobi(m)
    }

    pub fn extract(self) -> Eigen {
        let v = self.v;
        let d = self.d;
        Eigen {
            eigenvalue: d,
            eigenvector: v
        }
    }

    /// Main Jacobi algorithm
    ///
    /// * Reference : Press, William H., and William T. Vetterling. *Numerical Recipes.* Cambridge: Cambridge Univ. Press, 2007.
    pub fn iter(&mut self) {
        let a = &mut self.a;
        let n = self.n;
        let v = &mut self.v;
        let mut b = a.diag();
        let d = &mut self.d;
        let mut z = vec![0f64; n];
        let mut h: f64;
        let mut _n_rot = self.n_rot;

        for i in 1 .. 51 {
            let mut sm = 0f64;
            for ip in 0 .. n-1 {
                for iq in ip+1 .. n {
                    sm += a[(ip, iq)].abs();
                }
            }
            if sm == 0f64 {
                eigsrt(d, v);
                return ();
            }
            let tresh = if i < 4 {
                0.2 * sm / (n.pow(2) as f64)
            } else {
                0f64
            };
            for ip in 0 .. n-1 {
                for iq in ip+1 .. n {
                    let g = 100f64 * a[(ip, iq)].abs();
                    if i > 4 && g <= EPSILON * d[ip].abs() && g <= EPSILON * d[iq].abs() {
                        a[(ip, iq)] = 0f64;
                    } else if a[(ip, iq)].abs() > tresh {
                        h = d[iq] - d[ip];
                        let t = if g <= EPSILON * h.abs() {
                            a[(ip, iq)] / h
                        } else {
                            let theta = 0.5 * h / a[(ip, iq)];
                            let mut temp = 1f64 / (theta.abs() + (1f64 + theta.powi(2)).sqrt());
                            if theta < 0f64 {
                                temp = -temp;
                            }
                            temp
                        };
                        let c = 1f64 / (1f64 + t.powi(2)).sqrt();
                        let s = t * c;
                        let tau = s / (1f64 + c);
                        h = t * a[(ip, iq)];
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        a[(ip, iq)] = 0f64;
                        for j in 0 .. ip {
                            rot(a, s, tau, j, ip, j, iq);
                        }
                        for j in ip+1 .. iq {
                            rot(a, s, tau, ip, j, j, iq);
                        }
                        for j in iq+1 .. n {
                            rot(a, s, tau, ip, j, iq, j);
                        }
                        for j in 0 .. n {
                            rot(v, s, tau, j, ip, j, iq);
                        }
                        _n_rot += 1;
                    }
                }
            }
            for ip in 0 .. n {
                b[ip] += z[ip];
                d[ip] = b[ip];
                z[ip] = 0f64;
            }
        }
        assert!(false, "Too many iterations in routine jacobi");
    }
}

fn rot(a: &mut Matrix, s: f64, tau: f64, i: usize, j: usize, k: usize, l: usize) {
    let g = a[(i, j)];
    let h = a[(k, l)];
    a[(i, j)] = g - s * (h + g * tau);
    a[(k, l)] = h + s * (g - h * tau);
}

/// Given eigenvalue & eigenvector, sorts thod eigenvalues into descending order
///
/// * Reference : Press, William H., and William T. Vetterling. *Numerical Recipes.* Cambridge: Cambridge Univ. Press, 2007.
fn eigsrt(d: &mut Vec<f64>, v: &mut Matrix) {
    let mut k: usize;
    let n = d.len();
    for i in 0 .. n-1 {
        k = i;
        let mut p = d[k];
        for j in i .. n {
            if d[j] >= p {
                k = j;
                p = d[k];
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for j in 0 .. n {
                p = v[(j, i)];
                v[(j, i)] = v[(j, k)];
                v[(j, k)] = p;
            }
        }
    }
}