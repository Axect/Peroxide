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
struct RootFinder<F: Fn(f64) -> f64> {
    init: State,
    curr: State, 
    f: Box<F>,
}

//impl RootFinder<F: Fn(f64) -> f64> {
//    fn condition_number(&self) -> f64 {
//        match self.curr {
//            P(p) => {
//                let z = p;
//                let y = self.f(z);
//                let dy = 
//            }
//        }
//    }
//}
