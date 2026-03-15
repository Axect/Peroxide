//! Demonstrates the new Jet<N> const-generic automatic differentiation.
use peroxide::fuga::*;

fn main() {
    // 1. First-order: f(x) = x^3 at x = 2
    let x = Jet::<1>::var(2.0);
    let y = x.powi(3);
    println!("f(x) = x^3 at x = 2");
    println!("  f(2)  = {}", y.value());  // 8
    println!("  f'(2) = {}", y.dx());     // 12

    // 2. Second-order: f(x) = sin(x) at x = pi/4
    let x = Jet::<2>::var(std::f64::consts::FRAC_PI_4);
    let y = x.sin();
    println!("\nf(x) = sin(x) at x = pi/4");
    println!("  f(pi/4)   = {:.6}", y.value());  //  0.707107
    println!("  f'(pi/4)  = {:.6}", y.dx());     //  0.707107  (cos(pi/4))
    println!("  f''(pi/4) = {:.6}", y.ddx());    // -0.707107  (-sin(pi/4))

    // 3. Higher-order: f(x) = exp(x) at x = 0
    //    All derivatives of exp are 1 at 0, so f^(k)(0) = 1 for every k.
    let x = Jet::<5>::var(0.0);
    let y = x.exp();
    println!("\nf(x) = exp(x) at x = 0");
    for k in 0..=5 {
        println!("  f^({k})(0) = {:.1}", y.derivative(k));
    }

    // 4. Composition: f(x) = exp(sin(x)) at x = 0
    //    f(0)   = exp(0) = 1
    //    f'(0)  = cos(0)*exp(sin(0)) = 1
    //    f''(0) = (cos^2(0) - sin(0)) * exp(sin(0)) = 1
    //    f'''(0): can be computed but requires the chain rule repeatedly
    let x = Jet::<3>::var(0.0);
    let y = x.sin().exp();
    println!("\nf(x) = exp(sin(x)) at x = 0");
    println!("  f(0)    = {:.6}", y.value());
    println!("  f'(0)   = {:.6}", y.dx());
    println!("  f''(0)  = {:.6}", y.ddx());
    println!("  f'''(0) = {:.6}", y.derivative(3));

    // 5. Type aliases: Dual = Jet<1> and HyperDual = Jet<2>
    let x: Dual = Dual::var(1.0);
    let y = x.ln();
    println!("\nUsing Dual alias: f(x) = ln(x) at x = 1");
    println!("  f(1)  = {}", y.value());  // 0
    println!("  f'(1) = {}", y.dx());     // 1

    let x: HyperDual = HyperDual::var(1.0);
    let y = x.ln();
    println!("\nUsing HyperDual alias: f(x) = ln(x) at x = 1");
    println!("  f(1)   = {}", y.value()); // 0
    println!("  f'(1)  = {}", y.dx());    // 1
    println!("  f''(1) = {}", y.ddx());   // -1

    // 6. Backward-compat constructors (AD0 / AD1 / AD2)
    //    AD is an alias for Jet<2>; AD1 and AD2 set derivatives directly.
    let a = AD1(3.0, 1.0); // f(x) = x, evaluated at 3, derivative = 1
    let b = AD2(3.0, 1.0, 0.0);
    println!("\nBackward-compat AD constructors:");
    println!("  AD1(3, 1).value() = {}", a.value());
    println!("  AD1(3, 1).dx()    = {}", a.dx());
    println!("  AD2(3, 1, 0).ddx() = {}", b.ddx()); // 0
}
