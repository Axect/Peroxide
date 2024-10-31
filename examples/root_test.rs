use anyhow::Result;
use peroxide::fuga::*;

fn main() -> Result<()> {
    let problem = Cubic;
    let bisect = BisectionMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let newton = NewtonMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let false_pos = FalsePositionMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let secant = SecantMethod {
        max_iter: 100,
        tol: 1e-6,
    };
    let result_bisect = bisect.find(&problem)?;
    let result_newton = newton.find(&problem)?;
    let result_false_pos = false_pos.find(&problem)?;
    let result_secant = secant.find(&problem)?;
    println!("{:?}", result_bisect);
    println!("{:?}", result_newton);
    println!("{:?}", result_false_pos);
    println!("{:?}", result_secant);

    Ok(())
}

struct Cubic;

impl Cubic {
    fn eval(&self, x: Pt<1>) -> Result<Pt<1>> {
        Ok([(x[0] - 1f64).powi(3)])
    }
}

impl RootFindingProblem<1, 1, (f64, f64)> for Cubic {
    fn function(&self, x: Pt<1>) -> Result<Pt<1>> {
        self.eval(x)
    }
    fn initial_guess(&self) -> (f64, f64) {
        (0.0, 2.0)
    }
}

impl RootFindingProblem<1, 1, f64> for Cubic {
    fn function(&self, x: Pt<1>) -> Result<Pt<1>> {
        self.eval(x)
    }
    fn initial_guess(&self) -> f64 {
        0.0
    }
    fn derivative(&self, x: Pt<1>) -> Result<Jaco<1, 1>> {
        Ok([[3.0 * (x[0] - 1f64).powi(2)]])
    }
}
