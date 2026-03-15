//! Higher-order automatic differentiation with Jet<N>.
//! Demonstrates computing up to 10th-order derivatives.
use peroxide::fuga::*;

fn main() {
    // Taylor series of sin(x) around x = 0.
    // The k-th normalized Taylor coefficient is c_k = f^(k)(0) / k!.
    // For sin: c_0=0, c_1=1, c_2=0, c_3=-1/6, c_4=0, c_5=1/120, ...
    println!("Taylor coefficients of sin(x) at x = 0:");
    println!("  (c_k = f^(k)(0) / k!, i.e., normalized coefficients)");
    let x = Jet::<10>::var(0.0);
    let s = x.sin();
    for k in 0..=10 {
        let c = s.taylor_coeff(k);
        if c.abs() > 1e-15 {
            println!("  c_{k} = {c:+.10}");
        } else {
            println!("  c_{k} =  0");
        }
    }
    // Expected pattern: 0, +1, 0, -1/6≈-0.1667, 0, +1/120≈+0.0083, 0, -1/5040, 0, +1/362880, 0

    // Verify: k-th derivative of x^5 at x = 1.
    // f(x)   = x^5      → f(1)    = 1
    // f'(x)  = 5x^4     → f'(1)   = 5
    // f''(x) = 20x^3    → f''(1)  = 20
    // f'''(x)= 60x^2    → f'''(1) = 60
    // f^(4)  = 120x     → f^(4)(1)= 120
    // f^(5)  = 120      → f^(5)(1)= 120
    println!("\nDerivatives of x^5 at x = 1:");
    let x = Jet::<5>::var(1.0);
    let y = x.powi(5);
    for k in 0..=5 {
        println!("  f^({k})(1) = {:.0}", y.derivative(k));
    }

    // Compute an 8th-order Taylor polynomial of exp(x) at x = 0, evaluate at x = 0.5.
    // The Taylor coefficients of exp at 0 are all 1/k!, so taylor_coeff(k) = 1/k! for all k.
    println!("\nTaylor approximation of exp(x) at x = 0.5:");
    let x0 = Jet::<8>::var(0.0);
    let e = x0.exp();
    let h = 0.5_f64;
    let mut approx = 0.0_f64;
    for k in 0..=8 {
        approx += e.taylor_coeff(k) * h.powi(k as i32);
    }
    println!("  T_8(0.5) = {approx:.15}");
    println!("  exp(0.5) = {:.15}", 0.5_f64.exp());
    println!("  error    = {:.2e}", (approx - 0.5_f64.exp()).abs());

    // Demonstrate that derivative(k) gives the raw factorial-scaled value.
    // For f(x) = exp(x) at x = 0: derivative(k) should equal k! * taylor_coeff(k) = k!/k! = 1.
    println!("\nRaw derivatives of exp(x) at x = 0 (all should be 1):");
    let x = Jet::<6>::var(0.0);
    let y = x.exp();
    for k in 0..=6 {
        println!("  f^({k})(0) = {:.1}", y.derivative(k));
    }
}
