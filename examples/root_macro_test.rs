#[macro_use]
extern crate peroxide;
use anyhow::Result;
use peroxide::fuga::*;

fn main() -> Result<()> {
    let root_bisect = bisection!(f, (0.0, 2.0), 100, 1e-6)?;
    let root_newton = newton!(f, 0.0, 100, 1e-6)?;
    let root_false_pos = false_position!(f, (0.0, 2.0), 100, 1e-6)?;
    let root_secant = secant!(f, (0.0, 2.0), 100, 1e-6)?;

    let closure_bisect = bisection!(|x: f64| (x - 1f64).powi(3), (0.0, 2.0), 100, 1e-6)?;
    let closure_false_pos = false_position!(|x: f64| (x - 1f64).powi(3), (0.0, 2.0), 100, 1e-6)?;
    let closure_secant = secant!(|x: f64| (x - 1f64).powi(3), (0.0, 2.0), 100, 1e-6)?;

    println!("root_bisect: {}", root_bisect);
    println!("root_newton: {}", root_newton);
    println!("root_false_pos: {}", root_false_pos);
    println!("root_secant: {}", root_secant);

    println!("closure_bisect: {}", closure_bisect);
    println!("closure_false_pos: {}", closure_false_pos);
    println!("closure_secant: {}", closure_secant);

    assert!(f(root_bisect).abs() < 1e-6);
    assert!(f(root_newton).abs() < 1e-6);
    assert!(f(root_false_pos).abs() < 1e-6);
    assert!(f(root_secant).abs() < 1e-6);

    assert!(f(closure_bisect).abs() < 1e-6);
    assert!(f(closure_false_pos).abs() < 1e-6);
    assert!(f(closure_secant).abs() < 1e-6);

    Ok(())
}

#[ad_function]
fn f(x: f64) -> f64 {
    (x - 1f64).powi(3)
}
