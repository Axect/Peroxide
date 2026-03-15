extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    // The Real trait is implemented for f64 and AD (= Jet<2>).
    // Using the generic function with both types:
    let x_f64 = 2f64;
    let x_ad1 = AD1(2f64, 1f64);   // Jet<2> with first derivative = 1, second = 0
    let x_ad2 = AD2(2f64, 1f64, 0f64); // Jet<2> with explicit second derivative

    println!("f(x) = x^2, evaluated with Real trait:");
    println!("  f64:  f(2) = {}", f(x_f64));
    println!("  AD1:  f(2) = {}, f'(2) = {}", f(x_ad1).value(), f(x_ad1).dx());
    println!("  AD2:  f(2) = {}, f'(2) = {}, f''(2) = {}",
        f(x_ad2).value(), f(x_ad2).dx(), f(x_ad2).ddx());

    // Direct Jet<N> usage (without the Real trait — more explicit):
    println!("\nf(x) = x^2 with explicit Jet<N> types:");
    let x1 = Jet::<1>::var(2.0);
    let y1 = x1.powi(2);
    println!("  Jet<1>: f(2) = {}, f'(2) = {}", y1.value(), y1.dx());

    let x2 = Jet::<2>::var(2.0);
    let y2 = x2.powi(2);
    println!("  Jet<2>: f(2) = {}, f'(2) = {}, f''(2) = {}",
        y2.value(), y2.dx(), y2.ddx());

    let x3 = Jet::<3>::var(2.0);
    let y3 = x3.powi(2);
    println!("  Jet<3>: f(2) = {}, f'(2) = {}, f''(2) = {}, f'''(2) = {}",
        y3.value(), y3.dx(), y3.ddx(), y3.derivative(3));
}

// Generic function over the Real trait (works with f64 and AD = Jet<2>).
fn f<T: Real>(x: T) -> T {
    x.powi(2)
}
