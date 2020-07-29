use std::collections::HashMap;
use crate::structure::ad::{AD, AD1, AD2};
//use std::marker::PhantomData;
//use RootState::{P, I};

pub trait RootFinder {
    type InitialState;
    type RootType;

    fn mut_update(&mut self);
    fn find_root(&mut self) -> Option<Self::RootType>;
    fn set_initial_condition(&mut self, init: Self::InitialState) -> &mut Self;
    fn set_times(&mut self) -> &mut Self;
    fn check_enough(&self) -> bool;
}

#[derive(Debug, Clone, Copy, Hash, PartialOrd, PartialEq, Eq)]
enum RootOptions {
    InitCondition,
    Times,
}

pub struct Bisection<F: Fn(f64) -> f64> {
    init: (f64, f64),
    curr: (f64, f64),
    func: Box<F>,
    times: usize,
    options: HashMap<RootOptions, bool>,
}

//#[derive(Debug, Copy, Clone)]
//pub enum RootState {
//    P(f64),
//    I(f64, f64),
//}
//
///// Structure for Root finding
/////
///// * **Caution**: `T` $\geq$ `AD1`
//#[derive(Debug, Clone)]
//pub struct RootFinder<T, F> 
//where T: AD, F: Fn(T) -> T {
//    init: RootState,
//    pub curr: RootState, 
//    method: RootFind,
//    f: Box<F>,
//    _marker: PhantomData<T>
//}
//
//impl<T, F> RootFinder<T, F> 
//where T: AD, F: Fn(T) -> T {
//    pub fn new(init: RootState, method: RootFind, f: F) -> Self {
//        RootFinder {
//            init,
//            curr: init,
//            method,
//            f: Box::new(f),
//            _marker: PhantomData
//        }
//    }
//
//    pub fn condition_number(&self) -> f64 {
//        match self.curr {
//            P(p) => {
//                let mut z = AD1::from(p);
//                z.d1 = 1f64;
//                let fz = (self.f)(z.into()).to_ad1();
//                p * fz.d1 / fz.d0
//            }
//            I(a, b) => {
//                let p = (a + b) / 2f64;
//                let mut z = AD1::from(p);
//                z.d1 = 1f64;
//                let fz = (self.f)(z.into()).to_ad1();
//                p * fz.d1 / fz.d0
//            }
//        }
//    }
//
//    // pub fn root_find(&self) -> Option<f64> {
//    //     match self.method {
//    //         RootFind::Bisection => {
//    //
//    //         }
//    //     }
//    // }
//}
