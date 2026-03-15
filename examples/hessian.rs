use peroxide::fuga::*;

fn main() {
    // Compute partial derivatives of the Rosenbrock function
    // f(x, y) = (1 - x)^2 + 5*(y - x^2)^2
    //
    // Analytical values at (1, 1):
    //   f(1,1)      =  0
    //   df/dx       =  0
    //   df/dy       =  0
    //   d²f/dx²     =  42   (= 2*(1) + 5*2*(1*(-2x))^2|_{x=1} + ... = 2 + 40 = 42)
    //   d²f/dy²     =  10   (= 5*2 = 10)

    // df/dx and d²f/dx² at (1, 1): seed x as variable, keep y constant.
    let x = Jet::<2>::var(1.0);
    let y = Jet::<2>::constant(1.0);
    let result_x = f(x, y);
    println!("At (1, 1):");
    println!("  f(1, 1)   = {}", result_x.value());  // 0
    println!("  df/dx     = {}", result_x.dx());      // 0
    println!("  d²f/dx²   = {}", result_x.ddx());     // 42

    // df/dy and d²f/dy² at (1, 1): seed y as variable, keep x constant.
    let x = Jet::<2>::constant(1.0);
    let y = Jet::<2>::var(1.0);
    let result_y = f(x, y);
    println!("  df/dy     = {}", result_y.dx());      // 0
    println!("  d²f/dy²   = {}", result_y.ddx());     // 10
}

fn f(x: Jet<2>, y: Jet<2>) -> Jet<2> {
    let one = Jet::<2>::constant(1.0);
    let five = Jet::<2>::constant(5.0);
    (one - x).powi(2) + five * (y - x.powi(2)).powi(2)
}
