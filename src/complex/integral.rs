use crate::complex::C64;
use crate::numerical::integral::{GKIntegrable, GLKIntegrable, NCIntegrable};
use crate::structure::polynomial::{lagrange_polynomial, Calculus, Polynomial};

// Newton Cotes Quadrature for Complex Functions of one Real Variable
impl NCIntegrable for C64 {
    type NodeY = (Vec<f64>, Vec<f64>);
    type NCPolynomial = (Polynomial, Polynomial);

    fn compute_node_y<F>(f: F, node_x: &[f64]) -> Self::NodeY
    where
        F: Fn(f64) -> Self,
    {
        node_x
            .iter()
            .map(|x| {
                let z = f(*x);
                (z.re, z.im)
            })
            .unzip()
    }

    fn compute_polynomial(node_x: &[f64], node_y: &Self::NodeY) -> Self::NCPolynomial {
        (
            lagrange_polynomial(node_x.to_vec(), node_y.0.to_vec()),
            lagrange_polynomial(node_x.to_vec(), node_y.1.to_vec()),
        )
    }

    fn integrate_polynomial(p: &Self::NCPolynomial) -> Self::NCPolynomial {
        (p.0.integral(), p.1.integral())
    }

    fn evaluate_polynomial(p: &Self::NCPolynomial, x: f64) -> Self {
        p.0.eval(x) + C64::I * p.1.eval(x)
    }
}

// Gauss Lagrange and Kronrod Quadrature for Complex Functions of one Real Variable
impl GLKIntegrable for C64 {
    const ZERO: Self = C64::ZERO;
}

// Gauss Kronrod Quadrature for Complex Functions of one Real Variable
impl GKIntegrable for C64 {
    fn is_finite(&self) -> bool {
        C64::is_finite(*self)
    }

    fn gk_norm(&self) -> f64 {
        self.norm()
    }
}
