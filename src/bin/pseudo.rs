extern crate peroxide;
use peroxide::*;

fn main() {
    let a = matrix(
        vec![
            3.,
            1f64.sin(),
            1f64.sin(),
            8.,
            -1250.,
            2.,
            -(1f64.exp()),
            -(1f64.exp()),
            20.,
        ],
        3,
        3,
        Row,
    );

    a.print();

    let pqlu = a.lu().unwrap();
    let p = pqlu.p;
    let q = pqlu.q;
    let l = pqlu.l;
    let u = pqlu.u;

    println!("");
    println!("{:?}", p);
    println!("{:?}", q);
    println!("");

    let mut lu = l * u;
    lu.print();
    println!("");

    for (idx1, idx2) in q.into_iter().rev() {
        lu = lu.swap(idx1, idx2, Col);
    }
    lu.print();
    println!("");

    for (idx1, idx2) in p.into_iter().rev() {
        lu = lu.swap(idx1, idx2, Row);
    }
    lu.print();

    println!("");

    a.inv().unwrap().print();
}
