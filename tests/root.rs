#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;
use peroxide::numerical::root::*;

#[test]
fn test_root_find() -> Result<(), RootError> {
    let init = RootState::I(0.1f64, 5f64);
    let mut a1 = RootFinder::<AD1>::new(init, RootFind::Bisection, f_exp)?;
    let mut b1 = RootFinder::<AD1>::new(init, RootFind::Bisection, f_ln)?;
    let mut c1 = RootFinder::<AD1>::new(init, RootFind::Bisection, f_sqrt)?;
    let mut a2 = RootFinder::<AD1>::new(init, RootFind::FalsePosition, f_exp)?;
    let mut b2 = RootFinder::<AD1>::new(init, RootFind::FalsePosition, f_ln)?;
    let mut c2 = RootFinder::<AD1>::new(init, RootFind::FalsePosition, f_sqrt)?;
    a2.set_times(10000);
    b2.set_times(10000);
    c2.set_times(10000);
    let x1 = a1.find_root()?;
    let x2 = b1.find_root()?;
    let x3 = c1.find_root()?;
    let y1 = a2.find_root()?;
    let y2 = b2.find_root()?;
    let y3 = c2.find_root()?;
    x1.print();
    x2.print();
    x3.print();
    y1.print();
    y2.print();
    y3.print();
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
