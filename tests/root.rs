// #[macro_use]
extern crate peroxide;
use peroxide::fuga::*;
// use peroxide::numerical::root::*;

#[test]
fn test_root_finder() -> Result<(), RootError> {
    let init = RootState::I(0.1f64, 5f64);
    let init_p = RootState::P(0.1f64);
    let init_s = RootState::I(0.1f64, 0.2f64);
    // Bisection
    let mut a1 = RootFinder::<_, AD1>::new(init, RootFind::Bisection, f_exp)?;
    let mut b1 = RootFinder::<_, AD1>::new(init, RootFind::Bisection, f_ln)?;
    let mut c1 = RootFinder::<_, AD1>::new(init, RootFind::Bisection, f_sqrt)?;
    a1.set_tol(1e-15);
    b1.set_tol(1e-15);
    c1.set_tol(1e-15);
    // FalsePosition
    let mut a2 = RootFinder::<_, AD1>::new(init, RootFind::FalsePosition, f_exp)?;
    let mut b2 = RootFinder::<_, AD1>::new(init, RootFind::FalsePosition, f_ln)?;
    let mut c2 = RootFinder::<_, AD1>::new(init, RootFind::FalsePosition, f_sqrt)?;
    a2.set_times(10000);
    b2.set_times(10000);
    c2.set_times(10000);
    // Secant
    let mut a3 = RootFinder::<_, AD1>::new(init_s, RootFind::Secant, f_exp)?;
    let mut b3 = RootFinder::<_, AD1>::new(init_s, RootFind::Secant, f_ln)?;
    let mut c3 = RootFinder::<_, AD1>::new(init_s, RootFind::Secant, f_sqrt)?;
    a3.set_times(10000);
    b3.set_times(10000);
    c3.set_times(10000);
    // Newton
    let mut a4 = RootFinder::<_, AD1>::new(init_p, RootFind::Newton, f_exp)?;
    let mut b4 = RootFinder::<_, AD1>::new(init_p, RootFind::Newton, f_ln)?;
    let mut c4 = RootFinder::<_, AD1>::new(init_p, RootFind::Newton, f_sqrt)?;
    a4.set_tol(1e-15);
    b4.set_tol(1e-15);
    c4.set_tol(1e-15);

    let x1 = a1.find_root()?;
    let x2 = b1.find_root()?;
    let x3 = c1.find_root()?;
    let y1 = a2.find_root()?;
    let y2 = b2.find_root()?;
    let y3 = c2.find_root()?;
    let z1 = a3.find_root()?;
    let z2 = b3.find_root()?;
    let z3 = c3.find_root()?;
    let w1 = a4.find_root()?;
    let w2 = b4.find_root()?;
    let w3 = c4.find_root()?;
    x1.print();
    x2.print();
    x3.print();
    y1.print();
    y2.print();
    y3.print();
    z1.print();
    z2.print();
    z3.print();
    w1.print();
    w2.print();
    w3.print();

    Ok(())
}

#[test]
fn test_bisection() -> Result<(), RootError> {
    let x1 = bisection(f_exp, (0f64, 5f64), 100, 1e-15)?;
    let x2 = bisection(f_ln, (0.1f64, 5f64), 100, 1e-15)?;
    let x3 = bisection(f_sqrt, (0f64, 5f64), 100, 1e-15)?;
    let x4 = bisection(f_sin, (0f64, 5f64), 100, 1e-15)?;

    println!("Bisection: ");
    x1.print();
    x2.print();
    x3.print();
    x4.print();

    Ok(())
}

#[test]
fn test_secant() -> Result<(), RootError> {
    let x1 = secant(f_exp, (0f64, 0.1f64), 200, 1e-15)?;
    let x2 = secant(f_ln, (0.1f64, 0.2f64), 200, 1e-15)?;
    let x3 = secant(f_sqrt, (0f64, 0.1f64), 200, 1e-15)?;
    let x4 = secant(f_sin, (2f64, 2.1f64), 200, 1e-15)?;

    println!("Secant: ");
    x1.print();
    x2.print();
    x3.print();
    x4.print();

    Ok(())
}

#[test]
fn test_false_position() -> Result<(), RootError> {
    println!("False Position: ");
    let x1 = false_position(f_exp, (0f64, 5f64), 1000, 1e-10)?;
    let x2 = false_position(f_ln, (0.1f64, 5f64), 1000, 1e-10)?;
    let x3 = false_position(f_sqrt, (0f64, 5f64), 1000, 1e-10)?;
    let x4 = false_position(f_sin, (1f64, 5f64), 1000, 1e-10)?;

    x1.print();
    x2.print();
    x3.print();
    x4.print();

    Ok(())
}

#[test]
fn test_newton() -> Result<(), RootError> {
    let x1 = newton(f_exp, 0.1f64, 100, 1e-15)?;
    let x2 = newton(f_ln, 0.1f64, 100, 1e-15)?;
    let x3 = newton(f_sqrt, 0.1f64, 100, 1e-15)?;
    let x4 = newton(f_sin, 2f64, 100, 1e-15)?;

    println!("Newton: ");
    x1.print();
    x2.print();
    x3.print();
    x4.print();

    Ok(())
}

fn f_exp<T: AD>(x: T) -> T {
    x.exp() - 2f64
}

fn f_ln<T: AD>(x: T) -> T {
    x.ln()
}

fn f_sqrt<T: AD>(x: T) -> T {
    x.sqrt() - 2f64
}

fn f_sin<T: AD>(x: T) -> T {
    x.sin()
}
