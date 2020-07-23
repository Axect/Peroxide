extern crate peroxide;
use peroxide::fuga::*;

macro_rules! rnd {
    ( $x:expr ) => {
        ($x * 100.0).round() / 100.0
    };
}

fn f(x: f64) -> f64 {
    1.0 / (1.0 + x * x)
}

#[test]
fn test_cubic_spline_initialization() {
    let mut vx = Vec::new();
    let mut vy = Vec::new();
    for i in 0..11 {
        let x = i as f64;
        vx.push(x);
        vy.push(f(x));
    }

    let spline = CubicSpline::from_nodes(vx, vy);

    for i in 0..11 {
        let x = i as f64;
        let y = spline.eval(x);

        assert_eq!(rnd!(y), rnd!(f(x)));
    }
}

#[test]
fn test_cubic_spline_extension() {
    let mut vx = Vec::new();
    let mut vy = Vec::new();
    for i in 0..11 {
        let x = (i - 10) as f64;
        vx.push(x);
        vy.push(f(x));
    }

    let mut spline = CubicSpline::from_nodes(vx, vy);

    vx = Vec::new();
    vy = Vec::new();

    for i in 11..21 {
        let x = (i - 10) as f64;
        vx.push(x);
        vy.push(f(x));
    }

    spline.extend_with_nodes(vx, vy);

    for i in 0..21 {
        let x = (i - 10) as f64;
        let y = spline.eval(x);

        assert_eq!(
            format!("{} = {}", x, rnd!(y)),
            format!("{} = {}", x, rnd!(f(x)))
        );
    }
}
