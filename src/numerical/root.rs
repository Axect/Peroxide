use std::marker::PhantomData;
use crate::structure::ad::{AD, AD1, AD2};
use RootState::{P, I};

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum RootFind {
    Bisection,
    Newton,
    FalsePosition,
    Halley
}

#[derive(Debug, Copy, Clone)]
pub enum RootState {
    P(f64),
    I(f64, f64),
}

/// Structure for Rootfinding
///
/// * **Caution**: `T` $\geq$ `AD1`
#[derive(Debug, Clone)]
pub struct RootFinder<T, F> 
where T: AD, F: Fn(T) -> T {
    init: RootState,
    pub curr: RootState, 
    method: RootFind,
    f: Box<F>,
    _marker: PhantomData<T>
}

impl<T, F> RootFinder<T, F> 
where T: AD, F: Fn(T) -> T {
    pub fn new(init: RootState, method: RootFind, f: F) -> Self {
        RootFinder {
            init,
            curr: init,
            method,
            f: Box::new(f),
            _marker: PhantomData
        }
    }

    pub fn condition_number(&self) -> f64 {
        match self.curr {
            P(p) => {
                let mut z = AD1::from(p);
                z.d1 = 1f64;
                let fz = (self.f)(z.into()).to_ad1();
                p * fz.d1 / fz.d0
            }
            I(a, b) => {
                let p = (a + b) / 2f64;
                let mut z = AD1::from(p);
                z.d1 = 1f64;
                let fz = (self.f)(z.into()).to_ad1();
                p * fz.d1 / fz.d0
            }
        }
    }
}
