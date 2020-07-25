use std::marker::PhantomData;
use crate::structure::ad::{AD, AD1, AD2};
use State::{P, I};

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum RootFind {
    Bisection,
    Newton,
    FalsePosition,
}

#[derive(Debug, Copy, Clone)]
enum State {
    P(f64),
    I(f64, f64),
}

#[derive(Debug, Clone)]
struct RootFinder<T, F> 
where T: AD, F: Fn(T) -> T {
    init: State,
    curr: State, 
    f: Box<F>,
    _marker: PhantomData<T>
}

impl<T, F> RootFinder<T, F> 
where T: AD, F: Fn(T) -> T {
    fn condition_number(&self) -> f64 {
        match self.curr {
            P(p) => {
                let z = AD1::from(p);
                let fz = (self.f)(z.into()).to_ad1();
                p * fz.d1 / fz.d0
            }
            I(a, b) => {
                let p = (a + b) / 2f64;
                let z = AD1::from(p);
                let fz = (self.f)(z.into()).to_ad1();
                p * fz.d1 / fz.d0
            }
        }
    }
}
