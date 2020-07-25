#[macro_use]
extern crate peroxide;
use peroxide::fuga::*;
use peroxide::numerical::root::*;

#[test]
fn test_condition_number() {
    let zs = c!(1,2,3,4,5);
    let mut rf1 = RootFinder::<AD1, _>::new(RootState::P(1f64), RootFind::Bisection, f_exp);
    let mut rf2 = RootFinder::<AD1, _>::new(RootState::P(1f64), RootFind::Bisection, f_ln);
    let mut rf3 = RootFinder::<AD1, _>::new(RootState::P(1f64), RootFind::Bisection, f_sqrt);
    for z in zs {
        rf1.curr = RootState::P(z);
        rf2.curr = RootState::P(z);
        rf3.curr = RootState::P(z);
        rf1.condition_number().print();
        rf2.condition_number().print();
        rf3.condition_number().print();
    }
}

fn f_exp<T: AD>(x: T) -> T {
    x.exp()
}

fn f_ln<T: AD>(x: T) -> T {
    x.ln()
}

fn f_sqrt<T: AD>(x: T) -> T {
    x.sqrt()
}
