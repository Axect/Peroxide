#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;
use anyhow::Result;

fn main() -> Result<()> {
    let root_bisect = bisection!(f, (0.0, 2.0), 100, 1e-6);
    let root_newton = newton!(f, 0.0, 100, 1e-6);
    let root_false_pos = false_position!(f, (0.0, 2.0), 100, 1e-6);
    let root_secant = secant!(f, (0.0, 2.0), 100, 1e-6);

    println!("root_bisect: {}", root_bisect);
    println!("root_newton: {}", root_newton);
    println!("root_false_pos: {}", root_false_pos);
    println!("root_secant: {}", root_secant);

    assert!(f(root_bisect).abs() < 1e-6);
    assert!(f(root_newton).abs() < 1e-6);
    assert!(f(root_false_pos).abs() < 1e-6);
    assert!(f(root_secant).abs() < 1e-6);

    Ok(())
}

#[ad_function]
fn f(x: f64) -> f64 {
    (x - 1f64).powi(3)
}
