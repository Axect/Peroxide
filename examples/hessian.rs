use peroxide::fuga::*;

fn main() {
    let x = AD2(1f64, 1f64, 0f64);
    let y = AD2(1f64, 1f64, 0f64);
    let x0 = AD2(1f64, 0f64, 0f64);
    let y0 = AD2(1f64, 0f64, 0f64);

    let dy = f(x, y);
    dy.print();

    let dy0 = f(x0, y0);
    dy0.print();
}

fn f(x: AD, y: AD) -> AD {
    (1f64 - x).powi(2) + 5f64 * (y - x.powi(2)).powi(2)
}
