extern crate peroxide;
use peroxide::fuga::*;

fn main() {
    let a = integrate(|x| x.sin(), (1f64, 2f64), GaussLegendre(15));
    let b = integrate(|x| x.sin(), (1f64, 2f64), GaussLegendre(29));
    let c = integrate(|x| x.sin(), (1f64, 2f64), G7K15(1e-16, 20));
    let d = integrate(|x| x.sin(), (1f64, 2f64), G10K21(1e-16, 20));
    let e = integrate(|x| x.sin(), (1f64, 2f64), G15K31(1e-16, 20));
    let f = integrate(|x| x.sin(), (1f64, 2f64), G20K41(1e-16, 20));
    let g = integrate(|x| x.sin(), (1f64, 2f64), G25K51(1e-16, 20));
    let h = integrate(|x| x.sin(), (1f64, 2f64), G30K61(1e-16, 20));
    a.print();
    b.print();
    c.print();
    d.print();
    e.print();
    f.print();
    g.print();
    h.print();
}
